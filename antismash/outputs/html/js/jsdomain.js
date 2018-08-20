var jsdomain = {
    version: '0.0.1',
    label_height: 30,
    text_height: 14,
    unique_id: 0
};

jsdomain.drawDomains = function(id, region, height, width) {
    var container = d3.select('#' + id);
    var single_orf_height = height + jsdomain.label_height;
    container.selectAll("svg").remove();
    var chart = container.append("svg")
        .attr("height", single_orf_height * region["orfs"].length + 10)
        .attr("width", width);

    max_orf_len = 0;
    for (i=0; i < region["orfs"].length; i++) {
        max_orf_len = Math.max(max_orf_len, region["orfs"][i].sequence.length);
    }

    for (i=0; i < region["orfs"].length; i++) {
        var orf = region["orfs"][i];
        var idx = jsdomain.unique_id++;
        var offset = height/10;
        var x = d3.scale.linear()
          .domain([1, max_orf_len])
          .range([0, width]);
        chart.append("text")
            .text(orf['id'])
            .attr("x", 5)
            .attr("y", (single_orf_height * i) + jsdomain.label_height - 5)
            .attr("class", "jsdomain-orflabel");

        var group = chart.append("g");
        group.append("line")
          .attr("x1", 0)
          .attr("y1", (single_orf_height * i) + jsdomain.label_height + (height/2))
          .attr("x2", x(orf.sequence.length))
          .attr("y2", (single_orf_height * i) + jsdomain.label_height + (height/2))
          .attr("class", "jsdomain-line");

        group.selectAll("rect")
            .data(orf['domains'])
        .enter().append("rect")
            .attr("x", function(d) { return x(d.start) })
            .attr("y", (single_orf_height * i) + jsdomain.label_height)
            .attr("rx", 17)
            .attr("ry", 17)
            .attr("width", function(d) { return x(d.end) - x(d.start) })
            .attr("height", single_orf_height - jsdomain.label_height)
            .attr("id", function(d, j) { return "details-orf-" + idx + "-" + j + "-domain"})
            .attr("class", "jsdomain-domain ")
            .attr("fill", function(d) { return jsdomain.get_fill_color(d['type']);  })
            .attr("stroke", function(d) { return jsdomain.get_stroke_color(d['type']);  })
            .attr("stroke-width", 2);

        group.selectAll("text")
            .data(orf['domains'])
        .enter().append("text")
            .text(function(d) { return jsdomain.get_label(d['type']); })
            .attr("x", function(d) { return x(d.start) + (x(d.end) - x(d.start))/2 - this.getComputedTextLength()/2 })
            .attr("y", (single_orf_height * i) + jsdomain.label_height * 1.25  + height/2)
            .attr("id", function(d, j) { return "details-orf-" + idx + "-" + j + "-text"})
            .attr("class", "jsdomain-text")
            .attr("font-size", jsdomain.text_height)
            .attr("font-weight", "bold");

        toolgroup = container.append("div").attr("id", "details-orf-" + idx);
        toolgroup.selectAll("div")
            .data(orf['domains'])
        .enter().append("div")
            .attr("class", "jsdomain-tooltip")
            .attr("id", function(d, j) { return "details-orf-" + idx + "-" + j + "-tooltip"})
            .html(function(d) { return jsdomain.generateTooltip(d, orf) } );
    }
    jsdomain.init();
}

jsdomain.get_stroke_color = function(type) {
    switch (type) {
        case "AMP-binding":
        case "AOX":
            return "rgb(87,22,128)"
        case "PCP":
        case "ACP":
            return "rgb(11,78,199)"
        case "Cglyc":
        case "CXglyc":
        case "Condensation_DCL":
        case "Condensation_LCL":
        case "Condensation_Starter":
        case "Condensation_Dual":
        case "Heterocyclization":
            return "rgb(59,59,140)"
        case "Epimerization":
            return "rgb(59,59,140)"
        case "NRPS-COM_Nterm":
        case "NRPS-COM_Cterm":
        case "PKS_Docking_Nterm":
        case "PKS_Docking_Cterm":
        case "Trans-AT_docking":
            return "rgb(71,71,159)"
        case "Thioesterase":
        case "TD":
            return "rgb(119,3,116)"
        case "PKS_KS":
            return "rgb(9,179,9)"
        case "PKS_AT":
            return "rgb(221,6,6)"
        case "PKS_KR":
            return "rgb(10,160,76)"
        case "PKS_DH":
        case "PKS_DH2":
        case "PKS_DHt":
            return "rgb(186,103,15)"
        case "PKS_ER":
            return "rgb(12,161,137)"
        case "Aminotran_1_2":
        case "Aminotran_3":
        case "Aminotran_4":
        case "Aminotran_5":
        case "Polyketide_cyc2":
        default:
            return "rgb(147,147,147)"
    }
}

jsdomain.get_fill_color = function(type) {
    switch (type) {
        case "AMP-binding":
        case "AOX":
            return "rgb(188,127,245)"
        case "PCP":
        case "ACP":
            return "rgb(129,190,247)"
        case "Cglyc":
        case "CXglyc":
        case "Condensation_DCL":
        case "Condensation_LCL":
        case "Condensation_Starter":
        case "Condensation_Dual":
        case "Heterocyclization":
            return "rgb(129,129,247)"
        case "Epimerization":
            return "rgb(129,129,247)"
        case "NRPS-COM_Nterm":
        case "NRPS-COM_Cterm":
        case "PKS_Docking_Nterm":
        case "PKS_Docking_Cterm":
        case "Trans-AT_docking":
            return "rgb(128,128,245)"
        case "Thioesterase":
        case "TD":
            return "rgb(245,196,242)"
        case "PKS_KS":
            return "rgb(129,247,129)"
        case "PKS_AT":
            return "rgb(247,129,129)"
        case "PKS_KR":
            return "rgb(128,246,128)"
        case "PKS_DH":
        case "PKS_DH2":
        case "PKS_DHt":
            return "rgb(247,190,129)"
        case "PKS_ER":
            return "rgb(129,247,243)"
        case "Aminotran_1_2":
        case "Aminotran_3":
        case "Aminotran_4":
        case "Aminotran_5":
        case "Polyketide_cyc2":
        default:
            return "rgb(218,218,218)"
    }
}

jsdomain.get_label = function(type) {
    switch (type) {
        case "AMP-binding":
        case "AOX":
            return "A"
        case "PCP":
        case "ACP":
        case "NRPS-COM_Nterm":
        case "NRPS-COM_Cterm":
        case "PKS_Docking_Nterm":
        case "PKS_Docking_Cterm":
        case "Trans-AT_docking":
        case "Aminotran_1_2":
        case "Aminotran_3":
        case "Aminotran_4":
        case "Aminotran_5":
        case "Polyketide_cyc2":
            return ""
        case "Cglyc":
        case "CXglyc":
        case "Condensation_DCL":
        case "Condensation_LCL":
        case "Condensation_Starter":
        case "Condensation_Dual":
        case "Heterocyclization":
            return "C"
        case "Epimerization":
            return "E"
        case "Thioesterase":
            return "TE"
        case "PKS_KS":
            return "KS"
        case "PKS_AT":
            return "AT"
        case "PKS_KR":
            return "KR"
        case "PKS_DH":
        case "PKS_DH2":
            return "DH"
        case "PKS_DHt":
            return "DHt"
        case "PKS_ER":
            return "ER"
        default:
            return type.split('_')[0];
    }
}

jsdomain.generateTooltip = function(domain, orf) {
    html = domain['type'] + "<br>";
    html += "Location: " + domain['start'] + "-" + domain['end'] + " AA<br>";
    if (domain['napdoslink'].length > 0) {
        html += "<a href=\"" + domain['napdoslink'] + "\" target=\"_blank\">" + "Analyze with NaPDoS</a><br>";
    }
    html += "<a href=\"" + domain['blastlink'] + "\" target=\"_blank\">" + "Run BlastP on this domain</a><br>";
    if (domain['predictions'].length > 0) {
        html += "<dl><dt>Substrate predictions:</dt>";
        for (var i = 0; i < domain['predictions'].length; i++ ) {
            html += "<dd>-" + domain['predictions'][i][0] + ": ";
            html += domain['predictions'][i][1] + "</dd>";
        }
        html += "</dl>";
    }
    if (domain['type'] == 'PKS_KR') {
        html += "<dl><dt>" + domain['kr_stereo'] + "</dt></dl>";
    }
    html += 'AA sequence: <a href="javascript:copyToClipboard(' + "'" + domain['sequence'] + "'" + ')">Copy to clipboard</a><br>';
    if (domain['dna_sequence']){
        html += 'Nucleotide sequence: <a href="javascript:copyToClipboard(' + "'" + domain['dna_sequence'] + "'" + ')">Copy to clipboard</a><br>';
    }

    return html;
}

jsdomain.tooltip_handler = function(ev) {
    var id = $(this).attr("id").replace("-domain", "-tooltip");
    var tooltip = $("#"+id);

    if (jsdomain.active_tooltip) {
        jsdomain.active_tooltip.hide();
    }
    jsdomain.active_tooltip = tooltip;

    if (tooltip.css("display") == 'none') {
        var offset = $(this).offset();
        tooltip.css("left", offset.left + 10);
        var this_parent = $(this).parent();
        var top_offset = this_parent.height()/(this_parent.children('line').length * 2);
        tooltip.css("top", offset.top + 30);
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

jsdomain.init = function() {
    $(".jsdomain-domain").click(jsdomain.tooltip_handler);
    $(".jsdomain-text").click(function() {
        var id = $(this).attr("id").replace("-text", "-domain");
        $("#"+id).click();
    });
    $(".jsdomain-textarea").click(function(event) {
        event.stopPropagation();
    });
};
