(function($) {
  if (!$('html').hasClass('font-checked')) {
        $('html').append('<span class="test-font"/>');
        $('.test-font').html('&#xe000;').css({
            fontFamily: 'testfont',
            opacity: 1
        });
  }
  // Read the cookie set by the font file
  // Source: http://www.sitepoint.com/how-to-deal-with-cookies-in-javascript/
  function readCookie(name) {
    var nameEQ = name + "=";
    var ca = document.cookie.split(';');
    for(var i=0;i < ca.length;i++) {
      var c = ca[i];
      while (c.charAt(0)==' ') c = c.substring(1,c.length);
      if (c.indexOf(nameEQ) == 0) return c.substring(nameEQ.length,c.length);
    }
    return null;
  }

  $(window).bind("load", function() {
    // find font cookie
    fontdownload = readCookie("FontDownloaded");
    if (!fontdownload) {
      $('html').addClass('no-font font-checked');
      fontFallback();
    } else {
      $('html').addClass('font font-checked').removeClass('no-font');
    }

    // raise fontFallback event
    function fontFallback() {
      $.event.trigger({
        type: "fontFallback"
      });
    }
  });

})(jQuery);
;/**/
