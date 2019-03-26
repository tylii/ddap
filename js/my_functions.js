
//show rows that have cliked status and hide others rows
$(document).ajaxStop(function(){  // excecute when the page is ready
    $(".btn-group .btn").click(function(){
        var inputValue = $(this).find("input").val();
        if(inputValue != 'all'){
            var target = $('#Table1 tr[data-status="' + inputValue + '"]');
            $("#Table1 > tbody > tr").not(target).hide();
            target.fadeIn();
        } else {
            $("#Table1 > tbody > tr").fadeIn();
        }
    });
    // Changing the class of status label to support Bootstrap 4
    var bs = $.fn.tooltip.Constructor.VERSION;
    var str = bs.split(".");
    if(str[0] == 4){
        $(".label").each(function(){
            var classStr = $(this).attr("class");
            var newClassStr = classStr.replace(/label/g, "badge");
            $(this).removeAttr("class").addClass(newClassStr);
        });
    }
});

// add a navigation bar to the page
$(function(){
    $("#nav-bar").load("https://tylii.github.io/ddap/html/nav.html");
  });




