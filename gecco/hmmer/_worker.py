import ctypes
import multiprocessing
import queue
import threading
import typing

import psutil

from pyhmmer.easel import Alphabet, DigitalSequence
from pyhmmer.plan7 import Pipeline, TopHits, HMM


class _HMMPipelineThread(threading.Thread):
    @staticmethod
    def _none_callback(hmm: HMM, total: int) -> None:
        pass

    def __init__(
        self,
        contigs: typing.Iterable[typing.Collection[DigitalSequence]],
        query_queue: "queue.Queue[typing.Optional[typing.Tuple[int, HMM]]]",
        query_count: multiprocessing.Value,  # type: ignore
        hits_queue: "queue.PriorityQueue[typing.Tuple[int, typing.List[TopHits]]]",
        kill_switch: threading.Event,
        callback: typing.Optional[typing.Callable[[HMM, int], None]],
        options: typing.Dict[str, typing.Any],
    ) -> None:
        super().__init__()
        self.options = options
        self.pipeline = Pipeline(alphabet=Alphabet.amino(), **options)
        self.contigs = contigs
        self.query_queue = query_queue
        self.query_count = query_count
        self.hits_queue = hits_queue
        self.callback = callback or self._none_callback
        self.kill_switch = kill_switch
        self.error: typing.Optional[BaseException] = None

    def run(self) -> None:
        while not self.kill_switch.is_set():
            args = self.query_queue.get()
            if args is None:
                self.query_queue.task_done()
                return
            else:
                index, query = args
            try:
                self.process(index, query)
                self.query_queue.task_done()
            except BaseException as exc:
                self.error = exc
                self.kill()
                return

    def kill(self) -> None:
        self.kill_switch.set()

    def process(self, index, query):
        for hits in self.search(query):
            self.hits_queue.put(hits)
        self.callback(query, self.query_count.value)  # type: ignore

    def search(self, query: HMM) -> typing.Iterator[TopHits]:
        for contig in self.contigs:
            yield self.pipeline.search_hmm(query, contig)
            self.pipeline.clear()


class _HMMSearch(object):

    def __init__(
        self,
        queries: typing.Iterable[HMM],
        contigs: typing.Collection[typing.Collection[DigitalSequence]],
        cpus: int = 0,
        callback: typing.Optional[typing.Callable[[HMM, int], None]] = None,
        **options,
    ) -> None:
        self.queries = queries
        self.contigs = contigs
        self.cpus = cpus
        self.callback = callback
        self.options = options

    def _new_thread(
        self,
        query_queue: "queue.Queue[typing.Optional[typing.Tuple[int, HMM]]]",
        query_count: multiprocessing.Value,  # type: ignore
        hits_queue: "queue.PriorityQueue[typing.Tuple[int, TopHits]]",
        kill_switch: threading.Event,
    ) -> _HMMPipelineThread:
        return _HMMPipelineThread(
            self.contigs,
            query_queue,
            query_count,
            hits_queue,
            kill_switch,
            self.callback,
            self.options,
        )

    def _single_threaded(self) -> typing.Iterator[TopHits]:
        # create the queues to pass the HMM objects around, as well as atomic
        # values that we use to synchronize the threads
        query_queue = queue.Queue()  # type: ignore
        query_count = multiprocessing.Value(ctypes.c_ulong)
        hits_queue = queue.Queue()  # type: ignore
        kill_switch = threading.Event()

        # create the thread (to recycle code)
        thread = self._new_thread(query_queue, query_count, hits_queue, kill_switch)

        # process each HMM independently and yield the result
        # immediately
        for index, query in enumerate(self.queries):
            query_count.value += 1
            thread.process(index, query)
            while not hits_queue.empty():
                yield hits_queue.get_nowait()

    def _multi_threaded(self) -> typing.Iterator[TopHits]:
        # create the queues to pass the HMM objects around, as well as atomic
        # values that we use to synchronize the threads
        hits_queue = queue.Queue()  # type: ignore
        query_queue = queue.Queue()  # type: ignore
        query_count = multiprocessing.Value(ctypes.c_ulong)
        kill_switch = threading.Event()

        # create and launch one pipeline thread per CPU
        threads = []
        for _ in range(self.cpus):
            thread = self._new_thread(query_queue, query_count, hits_queue, kill_switch)
            thread.start()
            threads.append(thread)

        # queue the HMMs passed as arguments
        for index, query in enumerate(self.queries):
            query_count.value += 1
            query_queue.put((index, query))

        # poison-pill the queue so that threads terminate when they
        # have consumed all the HMMs
        for _ in threads:
            query_queue.put(None)

        # wait for all threads to be completed, and yield
        # results as soon as available otherwise
        while any(thread.is_alive() for thread in threads):
            while not hits_queue.empty():
                yield hits_queue.get_nowait()
            alive = next((thread for thread in threads if thread.is_alive()), None)
            if alive is not None:
                alive.join(timeout=0.1)

        # make sure threads didn't exit with an error, raise it otherwise
        for thread in threads:
            thread.join()
            if thread.error is not None:
                raise thread.error

        # yield remaining results
        while not hits_queue.empty():
            yield hits_queue.get_nowait()#[1]

    def run(self) -> typing.Iterator[TopHits]:
        if self.cpus == 1:
            return self._single_threaded()
        else:
            return self._multi_threaded()


def hmmsearch(
    queries: typing.Iterable[HMM],
    contigs: typing.Collection[typing.Collection[DigitalSequence]],
    cpus: int = 0,
    callback: typing.Optional[typing.Callable[[HMM, int], None]] = None,
    **options,  # type: typing.Any
) -> typing.Iterator[TopHits]:
    """Search HMM profiles against a sequence database.

    This function is a patched version from `pyhmmer.hmmer.hmmsearch` that
    enables processing several databases in parallel with the same query.

    """
    # count the number of CPUs to use
    _cpus = cpus if cpus > 0 else psutil.cpu_count(logical=False) or multiprocessing.cpu_count()
    runner = _HMMSearch(queries, contigs, _cpus, callback, **options)
    return runner.run()
