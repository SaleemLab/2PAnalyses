function removeLastSeven(inputFile, outputFile)
    % Import the Java PDFBox library bundled with MATLAB
    import org.apache.pdfbox.pdmodel.PDDocument;
    import java.io.File;

    % Load the existing PDF
    pdfDoc = PDDocument.load(File(inputFile));
    
    % Get the current total number of pages
    totalRef = pdfDoc.getNumberOfPages();
    
    % Safety: Ensure we don't try to remove more pages than exist
    numToRemove = min(7, totalRef);
    
    try
        % Loop exactly 7 times (or totalRef times if the doc is tiny)
        for i = 1:numToRemove
            % Always remove the current last page
            pdfDoc.removePage(pdfDoc.getNumberOfPages() - 1);
        end
        
        % Save and free up the file lock
        pdfDoc.save(outputFile);
        pdfDoc.close();
        fprintf('Done. Removed the last %d pages from the document.\n', numToRemove);
        
    catch ME
        pdfDoc.close(); 
        rethrow(ME);
    end
end