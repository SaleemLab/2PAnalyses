function removeSecondHalf(inputFile, outputFile)
    import org.apache.pdfbox.pdmodel.PDDocument;
    import java.io.File;

    % Load the document
    pdfDoc = PDDocument.load(File(inputFile));
    
    % Get total page count
    totalRef = pdfDoc.getNumberOfPages();
    
    % Calculate starting point for removal (the second half)
    % Using floor ensures we keep the middle page if the count is odd
    startPage = floor(totalRef / 2) + 1; 
    
    try
        % We iterate backwards from the last page to the midpoint
        for p = totalRef:-1:startPage
            % PDFBox is 0-indexed, so page 10 is index 9
            pdfDoc.removePage(p - 1);
        end
        
        pdfDoc.save(outputFile);
        pdfDoc.close();
        fprintf('Done! Kept the first %d pages.\n', startPage - 1);
        
    catch ME
        pdfDoc.close();
        rethrow(ME);
    end
end