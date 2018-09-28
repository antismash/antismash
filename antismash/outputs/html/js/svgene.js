/* Copyright 2012 Kai Blin. Licensed under the Apache License v2.0, see LICENSE file */

var svgene = {
    version: "0.1.6",
    label_height: 14,
    extra_label_width: 100,
    unique_id: 0
};

svgene.geneArrowPoints = function (orf, height, offset, border, scale) {
  var top_ = offset + svgene.label_height + border;
  var bottom = offset + svgene.label_height + height - border;
  var middle = offset + svgene.label_height + (height/2);
  if (orf.strand == 1) {
      var start = scale(orf.start);
      var box_end = Math.max(scale(orf.end) - (2*border), start);
      var point_end = scale(orf.end);
      points  = "" + start + "," + top_;
      points += " " + box_end + "," + top_;
      points += " " + point_end + "," + middle;
      points += " " + box_end + "," + bottom;
      points += " " + start + "," + bottom;
      return points;
  }
  if (orf.strand == -1) {
      var point_start = scale(orf.end);
      var end = scale(orf.start);
      var box_start = Math.min(scale(orf.end) + (2*border), end);
      points = "" + point_start + "," + middle;
      points += " " + box_start + "," + top_;
      points += " " + end + "," + top_;
      points += " " + end + "," + bottom;
      points += " " + box_start + "," + bottom;
      return points;
  }
};

svgene.ttaCodonPoints = function (codon, height, offset, border, scale) {
    var top_ = offset + svgene.label_height + height;;
    var bottom = offset + (2 * svgene.label_height) + height - border;
    var tip = Math.floor(scale(codon.start), scale(codon.end));
    var points = "" + tip + "," + top_;
    points += " " + (tip - 5) + "," + bottom;
    points += " " + (tip + 5) + "," + bottom;
    return points;
}

svgene.drawOrderedRegionOrfs = function(region, chart, all_orfs, borders, tta_codons,
                                         scale, idx, height, width, offset) {
  /* region: the region JSON object created in python
     chart:
     all_orfs: the ORF JSON objects created in python
     borders: the ClusterBorder JSON objects created in python
     tta_codons: the TTA codon JSON objects created in python
     scale:
     idx:
     height: the draw height for ORFs
     width: the draw width of the ORF centerline
     offset:
  */
  var orf_y = borders.length * 12 + svgene.label_height;

  // ORF centerline
  chart.append("line")
    .attr("x1", 0)
    .attr("y1", orf_y + svgene.label_height + height / 2)
    .attr("x2", width)
    .attr("y2", orf_y + svgene.label_height + height / 2)
    .attr("class", "svgene-line");

  // clusters
  var bar_size = 10;
  var vertical_bar_gap = 2;
  var cluster_bars = chart.selectAll("g")
    .data(borders)
  .enter().append("g")
    .attr("transform", function(d, i) { return "translate(0,0)"});
  // extent lines
  cluster_bars.append("line")
    .attr("x1", function(d){ return scale(d.neighbouring_start); })
    .attr("y1", function(d){ return d.height * (bar_size + vertical_bar_gap) + offset + bar_size/2})
    .attr("x2", function(d){ return scale(d.neighbouring_end); })
    .attr("y2", function(d){ return d.height * (bar_size + vertical_bar_gap) + offset + bar_size/2})
    .attr("class", "svgene-line");
  // rect containing first and last ORF triggering border detection
  cluster_bars.append("rect")
    .attr("width", function(d){ return Math.abs(scale(d.end) - scale(d.start))})
    .attr("height", bar_size)
    .attr("x", function(d){ return Math.min(scale(d.start), scale(d.end))})
    .attr("y", function(d){ return d.height * (bar_size + vertical_bar_gap) + offset})
    .attr("class", function(d){ return "svgene-border-" + d.tool});
  // cluster name
  cluster_bars.append("text")
    .attr("x", function(d) { return Math.min(scale(d.start), scale(d.end)) + 1 })
    .attr("y", function(d) { return (d.height + 2) * (bar_size) + d.height + offset - (2 - d.height)})
    .attr("dy", "-1em")
    .style("font-size", "xx-small")
    .attr("class", "clusterlabel")
    .text(function(d) { return d.product; });

  // ORFs
  chart.selectAll("polygon")
    .data(all_orfs)
  .enter().append("polygon")
    .attr("points", function(d) { return svgene.geneArrowPoints(d, height, orf_y, offset, scale); })
    .attr("class", function(d) { return "svgene-type-" + d.type + " svgene-orf"; })
    .attr("id", function(d) { return idx + "-region" + region.idx + "-" + svgene.tag_to_id(d.locus_tag) + "-svgeneorf"; })
    .attr("style", function(d) { if (d.color !== undefined) { return "fill:" + d.color; } });

  // TTA codons
  chart.selectAll("polyline.svgene-tta-codon")
    .data(tta_codons)
  .enter().append("polyline")
    .attr("points", function(d) { return svgene.ttaCodonPoints(d, height, orf_y, offset, scale) })
    .attr("class", "svgene-tta-codon");

  // ORF labels
  chart.selectAll("text.svgene-locustag")
    .data(all_orfs)
  .enter().append("text")
    .attr("x", function(d) { return scale(d.start); })
    .attr("y", orf_y + svgene.label_height)
    .attr("class", "svgene-locustag")
    .attr("id", function(d) { return idx + "-region" + region.idx + "-" + svgene.tag_to_id(d.locus_tag) + "-label"; })
    .text(function(d) { return d.locus_tag; });
};

svgene.drawRegion = function(id, region, height, width) {
  var container = d3.select("#" + id);
  container.selectAll("svg").remove();
  container.selectAll("div").remove();
  var chart = container.append("svg")
    .attr("height", 2 * height + (2 * svgene.label_height) + region.clusters.length * 12)
    .attr("width", width + svgene.extra_label_width);

  var all_orfs = [];
  var all_borders = [];
  var all_ttas = [];
  all_orfs.push.apply(all_orfs, region.orfs.sort(svgene.sort_biosynthetic_orfs_last));
  all_borders.push.apply(all_borders, region.clusters ? region.clusters : []);
  all_ttas.push.apply(all_ttas, region.tta_codons ? region.tta_codons: []);

  var idx = svgene.unique_id++;
  var offset = height/10;
  var x = d3.scale.linear()
    .domain([region.start, region.end])
    .range([0, width]);
  // only in rare cases, reassign the scale to a divergent function
  if (region.start > region.end) {
      x = function(position) {
          // pre-ori scale
          if (position < region.start) {
            return d3.scale.linear()
                  .domain([1, region.end])
                  .range([(region.sequence_length - region.start + region.end)/width, width])(position);
          }
          // post-ori scale
          else {
            return d3.scale.linear()
                  .domain([region.start, region.sequence_length])
                  .range([0, (region.sequence_length - region.start + region.end)/width])(position);
          }
      };
  }
  svgene.drawOrderedRegionOrfs(region, chart, all_orfs, all_borders, all_ttas,
                                x, idx, height, width, offset);
  container.selectAll("div")
    .data(all_orfs)
  .enter().append("div")
    .attr("class", "svgene-tooltip")
    .attr("id", function(d) { return idx + "-region" + region.idx + "-" + svgene.tag_to_id(d.locus_tag) + "-tooltip"; })
    .html(function(d) { return d.description});

  if (region.label !== undefined) {
    chart.append("text")
        .text(region.label)
        .attr("class", "svgene-regionlabel")
        .attr("x", function() { return width + svgene.extra_label_width - this.getComputedTextLength() - 5})
        .attr("y", function() { return svgene.label_height } )
        .attr("font-size", svgene.label_height);
  }

  svgene.init();
};

svgene.sort_biosynthetic_orfs_last = function(a, b) {
    if ((a.type != "biosynthetic" && b.type != "biosynthetic") ||
        (a.type == "biosynthetic" && b.type == "biosynthetic")) {
        return a.start - b.start;
    };
    if (a.type == "biosynthetic") {
        return 1;
    }
    return -1;
};

svgene.tag_to_id = function(tag) {
    return tag.replace(/(:|\.)/g, '-').replace(/-svgeneorf/g, '_orf');
}


svgene.tooltip_handler = function(ev) {
    var id = $(this).attr("id").replace("-svgeneorf", "-tooltip");
    var tooltip = $("#"+id);

    if (svgene.active_tooltip) {
        svgene.active_tooltip.hide();
    }
    svgene.active_tooltip = tooltip;

    if (tooltip.css("display") == 'none') {
        var offset = $(this).offset();
        tooltip.css("left", offset.left + 10);
        var this_parent = $(this).parent();
        var top_offset = this_parent.height()/(this_parent.children('line').length * 2);
        tooltip.css("top", offset.top + top_offset);
        tooltip.show();
        tooltip.click(function(){$(this).hide()});
        var timeout = setTimeout(function(){ tooltip.slideUp("fast") }, 5000);
        tooltip.data("timeout", timeout);
        tooltip.mouseover(function() {
            clearTimeout(tooltip.data("timeout"));
        }).mouseout(function() {
            timeout = setTimeout(function(){ tooltip.slideUp("fast") }, 5000);
            tooltip.data("timeout", timeout);
        });
    } else {
        tooltip.hide();
    }
};

svgene.init = function() {
    $(".svgene-orf").mouseover(function(e) {
        var id = $(this).attr("id").replace("-svgeneorf", "-label");
        $("#"+id).show();
    }).mouseout(function(e) {
        var id = $(this).attr("id").replace("-svgeneorf", "-label");
        $("#"+id).hide();
    }).click(svgene.tooltip_handler);
    $(".svgene-textarea").click(function(event) {
        event.stopPropagation();
    });
};
