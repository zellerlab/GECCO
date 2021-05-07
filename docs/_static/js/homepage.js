$(document).ready(function() {
  (function ($) {
    if (window.location.href.match("/index.html") !== null) {
      $("#gecco").children("p").last().html( (index, text) => text.replaceAll("</a> ", "</a>") ).end();
    }
  })(window.$jqTheme || window.jQuery);
})
