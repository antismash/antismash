/**
 * Created by Xiaowen on 12/23/16.
 * Adapted from clusterblast.js
 */

var domainalign = {
};

domainalign.init = function(parent_id) {
    $('#' + parent_id + ' .domainalign-orf').each(function() {
        var orf = $(this);
        domainalign.setLabel(orf, parent_id);
        domainalign.setTooltip(orf, parent_id);
    });
}

domainalign.setLabel = function(orf, parent_id) {
    var id =  orf.attr('id');
    var label = $('<div>');
    label.addClass('domainalign-locustag');
    label.attr('id', id + '-label');
    label.html(orf.attr('locus_tag').replace('[br]', '<br>'));
    // label.text(orf.attr('locus_tag'));

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

domainalign.setTooltip = function(orf, parent_id) {
    var id =  orf.attr('id');

    var tooltip = $('<div>');
    tooltip.addClass('domainalign-tooltip');
    tooltip.attr('id', id + '-tooltip');
    tooltip.html(orf.attr('description').replace(/\[br\]/g, '<br>'));
    $('#' + parent_id).append(tooltip);
    orf.click(domainalign.tooltip_handler);
}

domainalign.tooltip_handler = function(ev) {
    var orf_id = $(this).attr("id");
    var id = orf_id + "-tooltip";
    var tooltip = $("#"+id);

    if (domainalign.active_tooltip) {
        domainalign.active_tooltip.hide();
    }
    domainalign.active_tooltip = tooltip;

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

