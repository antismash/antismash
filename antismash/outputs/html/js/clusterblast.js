var clusterblast = {
};

clusterblast.init = function(parent_id) {
    $('#' + parent_id + ' .clusterblast-orf').each(function() {
        var orf = $(this);
        clusterblast.setLabel(orf, parent_id);
        clusterblast.setTooltip(orf, parent_id);
    });
}

clusterblast.setLabel = function(orf, parent_id) {
    //parent_id = clusterblast-16-svg
    var id =  orf.attr('id');
    var label = $('<div>');
    label.addClass('clusterblast-locustag');
    label.attr('id', id + '-label');
    label.text(orf.attr('locus_tag'));

    $('#' + parent_id).append(label);

    orf.mouseover(function(e) {
        var ofs = orf.offset();
        label.css('top', ofs.top - 32);
        label.css('left', ofs.left);
        $("#"+id+'-label').show();
    }).mouseout(function(e) {
        $("#"+id+'-label').hide();
    })
}

clusterblast.setTooltip = function(orf, parent_id) {
    var id =  orf.attr('id');

    var tooltip = $('<div>');
    tooltip.addClass('clusterblast-tooltip');
    tooltip.attr('id', id + '-tooltip');
    tooltip.html(orf.attr('description').replace('[br]', '<br>'));
    $('#' + parent_id).append(tooltip);
    orf.click(clusterblast.tooltip_handler);
}

clusterblast.tooltip_handler = function(ev) {
    var orf_id = $(this).attr("id");
    var id = orf_id + "-tooltip";
    var tooltip = $("#"+id);

    if (clusterblast.active_tooltip) {
        clusterblast.active_tooltip.hide();
    }
    clusterblast.active_tooltip = tooltip;

    if (tooltip.css("display") == 'none') {
        var ofs = $(this).offset();
        tooltip.css('top', ofs.top + 10);
        tooltip.css('left', ofs.left + 5);

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
