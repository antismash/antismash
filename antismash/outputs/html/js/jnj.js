
function toggle_downloadmenu(event) {
  event.preventDefault();
  $("#downloadmenu").fadeToggle("fast", "linear");
}

function switch_to_region() {
  setTimeout(function() {
    var url = $.url();
    $(".page").hide();
    $("li.regbutton").removeClass("active");
    var anchor = url.data.attr.fragment;
    if (anchor == "") {
      anchor = "overview";
    }
    $("#" + anchor).show();
    if (anchor != "overview") {
      $("li.regbutton." + anchor).addClass("active");
    }

    if (all_regions[anchor] !== undefined) {
      svgene.drawRegion(anchor+"-svg", all_regions[anchor], 20, 700);
    }
    if ($("#" + anchor + "-details-svg").length > 0) {
      jsdomain.drawDomains(anchor+ "-details-svg", details_data[anchor], 40, 700);
    }
    $("#" + anchor + " .clusterblast-selector").change();
  }, 1);
}

function next_region() {
  var regions = all_regions['order'];
  var current = $.url().data.attr.fragment;
  if (current == "" || current == "overview") {
    next = "r0c1";
  } else {
    current_index = regions.indexOf(current);
    if (current_index == regions.length - 1) {
      next = "overview";
    } else {
      next = regions[current_index + 1];
    }
  }
  window.location.href = "#" + next;
  switch_to_region();
}

function previous_region() {
  var regions = all_regions['order'];
  var current = $.url().data.attr.fragment;
  if (current == "" || current == "overview") {
    prev = regions[regions.length - 1];
  } else {
    current_index = regions.indexOf(current);
    if (current_index == 0) {
      prev = "overview";
    } else {
      prev = regions[current_index - 1];
    }
  }
  window.location.href = "#" + prev;
  switch_to_region();
}

function keyUpEvent(event) {
    var key = event.keyCode;
    if (key == 37) {  // left arrow
        previous_region();
    } else if (key == 39) {  // right arrow
        next_region();
    }
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
  document.addEventListener('keyup', keyUpEvent, false);
  $("#download").click(toggle_downloadmenu);

  $("#next-region").click(next_region);
  $("#prev-region").click(previous_region);

  $(".regbutton").click(function() {
    /* Make sure that even if user missed the link and clicked the
    background we still have the correct anchor */
    var href = $(this).children().first().attr('href');

    if (href === undefined) {
      return;
    }
    window.location.href = href;

    switch_to_region();
  }).mouseover(function() {
    /* Set the select region label text to region type */
    var classes = $(this).attr('class').split(' ');
    if (classes.length < 2) {
      return;
    }
    if (classes[1] == 'separator') {
      return;
    }
    var region_type = map_type_to_desc(classes[1]);
    var label = $('#region-type');
    label.data("orig_text", label.text());
    label.text(region_type + ":");
  }).mouseout(function() {
    /* and reset the select region label text */
    var label = $('#region-type');
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

  switch_to_region();
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
