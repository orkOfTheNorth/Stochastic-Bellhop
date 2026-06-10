function saveFig(fig, base_path)
% Save figure as 300-dpi PNG and (if available) PDF.
    print(fig, base_path, '-dpng', '-r300');
    try
        exportgraphics(fig, [base_path '.pdf'], 'ContentType', 'image', 'Resolution', 300);
    catch; end
end
