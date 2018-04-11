
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
  var clusters = geneclusters['order'];
  var current = $.url().data.attr.fragment;
  if (current == "" || current == "overview") {
    next = "r0c1";
  } else {
    current_index = clusters.indexOf(current);
    if (current_index == clusters.length - 1) {
      next = "overview";
    } else {
      next = clusters[current_index + 1];
    }
  }
  window.location.href = "#" + next;
  switch_to_cluster();
}

function previous_cluster() {
  var clusters = geneclusters['order'];
  var current = $.url().data.attr.fragment;
  if (current == "" || current == "overview") {
    prev = clusters[clusters.length - 1];
  } else {
    current_index = clusters.indexOf(current);
    if (current_index == 0) {
      prev = "overview";
    } else {
      prev = clusters[current_index - 1];
    }
  }
  window.location.href = "#" + prev;
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
  draw_structures();
});

function draw_structures () {

  $('.smiles-canvas').each(  function(counter) {
    //Draw structure on a 250x250 canvas
    let options = {
        width: 250,
        height: 250,
    };
    let smilesDrawer = new SmilesDrawer.Drawer(options);
    let canvas = this;
    SmilesDrawer.parse(canvas.getAttribute('data-smiles'), function(tree) {
            smilesDrawer.draw(tree, canvas, 'light', false);
        });
    //Calculate extra vertical space left for this structure
    drawing_height = smilesDrawer.canvasWrapper.drawingHeight;
    drawing_width = smilesDrawer.canvasWrapper.drawingWidth;
    vertical_size = (drawing_height/drawing_width) * 250;
    options = {
        width: 250,
        height: vertical_size * 0.8 //remove a bit of extra white space
    };
    //Draw it again; THIS IS NOT VERY ELEGANT, BUT WELL...
    smilesDrawer = new SmilesDrawer.Drawer(options);
    canvas.height = vertical_size;
    SmilesDrawer.parse(canvas.getAttribute('data-smiles'), function(tree) {
            smilesDrawer.draw(tree, canvas, 'light', false);
        });
    })
  
}
