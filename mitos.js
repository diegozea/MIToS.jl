$(document).ready(function(){

  // http://stackoverflow.com/questions/17341122/link-and-execute-external-javascript-file-hosted-on-github
  var raw = "https://cdn.rawgit.com/diegozea/MIToS.jl/master/docs/";

  var navcol_source   = $("#navcol-template").html();
  var navcol_template = Handlebars.compile(navcol_source);

  var content_source   = $("#content-template").html();
  var content_template = Handlebars.compile(content_source);

  var navcol_html     = navcol_template(pages);
  $("#NavCol").html(navcol_html)

  fillcontent = function(direction){
        $("#includedContent").load(direction);
  };

  fillcontent("home.html");

  $(".link").click(function(){

    $(".link").removeClass("active");
    $(this).addClass("active");

    var direction = $(this).data("src");

    fillcontent(direction);
  });

});
