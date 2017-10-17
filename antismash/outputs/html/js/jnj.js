
function toggle_downloadmenu(event) {
  event.preventDefault();
  $("#downloadmenu").fadeToggle("fast", "linear");
}

function switch_to_cluster() {
  setTimeout(function() {
    var url = $.url();
    $(".page").hide();
    $("li.clbutton").removeClass("active");
    var anchor = url.data.attr.fragment;
    if (anchor == "") {
      anchor = "overview";
    }
    $("#" + anchor).show();
    if (anchor != "overview") {
      $("li.clbutton." + anchor).addClass("active");
    }

    if (geneclusters[anchor] !== undefined) {
      svgene.drawClusters(anchor+"-svg", [geneclusters[anchor]], 20, 700);
    }
    if ($("#" + anchor + "-details-svg").length > 0) {
      jsdomain.drawDomains(anchor+ "-details-svg", details_data[anchor], 40, 700);
    }
    $("#" + anchor + " .clusterblast-selector").change();
  }, 1);
}

function next_cluster() {
  var num_clusters = Object.keys(geneclusters).length;
  var url = $.url();
  var anchor = url.data.attr.fragment;
  var href = "#" + anchor;
  if (anchor == "" || anchor == "overview") {
    anchor = "cluster-0";
  }
  var cluster_number = parseInt(anchor.split('-')[1]);
  var next_cluster_number = cluster_number + 1;
  if (next_cluster_number <= num_clusters) {
    href = "#cluster-" + next_cluster_number;
  } else {
    href = "#overview";
  }
  window.location.href = href;
  switch_to_cluster();
}

function previous_cluster() {
  var num_clusters = Object.keys(geneclusters).length;
  var url = $.url();
  var anchor = url.data.attr.fragment;
  var href = "#" + anchor;
  if (anchor == "" || anchor == "overview") {
    anchor = "cluster-0";
  }
  var cluster_number = parseInt(anchor.split('-')[1]);
  var prev_cluster_number = cluster_number - 1;
  if (prev_cluster_number == 0 ) {
    href = "#overview";
  } else if (prev_cluster_number < 0){
    href = "#cluster-" + num_clusters;
  } else {
    href = "#cluster-" + prev_cluster_number;
  }
  window.location.href = href;
  switch_to_cluster();
}

function toggle_cluster_rules(ev) {
  ev.preventDefault();
  var id = $(this).attr('id').replace(/-header/, '');
  var rules = $('#' + id);
  if (rules.css('display') == "none") {
    $(this).text('Hide pHMM detection rules used');
  } else {
    $(this).text('Show pHMM detection rules used');
  }
  rules.fadeToggle("fast", "linear");
}

function map_type_to_desc(type) {
  switch(type) {
    case "nrps": return "NRPS";
    case "t1pks": return "Type I PKS";
    case "t2pks": return "Type II PKS";
    case "t3pks": return "Type III PKS";
    case "t4pks": return "Type IV PKS";
    default: return type;
  }
}

function copyToClipboard (text) {
  window.prompt ("Copy to clipboard: Ctrl+C, Enter", text);
}

$(document).ready(function() {

  $("#download").click(toggle_downloadmenu);

  $("#next-cluster").click(next_cluster);
  $("#prev-cluster").click(previous_cluster);

  $(".clbutton").click(function() {
    /* Make sure that even if user missed the link and clicked the
    background we still have the correct anchor */
    var href = $(this).children().first().attr('href');

    if (href === undefined) {
      return;
    }
    window.location.href = href;

    switch_to_cluster();
  }).mouseover(function() {
    /* Set the select cluster label text to cluster type */
    var classes = $(this).attr('class').split(' ');
    if (classes.length < 2) {
      return;
    }
    if (classes[1] == 'separator') {
      return;
    }
    var cluster_type = map_type_to_desc(classes[1]);
    var label = $('#cluster-type');
    label.data("orig_text", label.text());
    label.text(cluster_type + ":");
  }).mouseout(function() {
    /* and reset the select cluster label text */
    var label = $('#cluster-type');
    label.text(label.data("orig_text"));
  });

  $('.clusterblast-selector').change(function() {
    var id = $(this).attr('id').replace('-select', '');
    var url = $(this).val();
    $.get(url, function(data) {
      $('#' + id + '-svg').html(data);
      clusterblast.init(id + '-svg');
      //            id =
    }, 'html');
    $('#' + id + '-download').off('click');
    $('#' + id + '-download').click(function () {
      var url = $("#" + id + "-select").val();
      window.open(url, '_blank');
    });
  });

  $('.cluster-rules-header').click(toggle_cluster_rules);

  switch_to_cluster();

});
