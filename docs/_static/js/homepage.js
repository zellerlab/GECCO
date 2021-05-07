$(document).ready(function() {
  (function ($) {
    $("#gecco").children("p").last().html( (index, text) => text.replaceAll("</a> ", "</a>") ).end();
  })(window.$jqTheme || window.jQuery);
})
