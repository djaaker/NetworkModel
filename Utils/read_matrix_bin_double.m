function matrix = read_matrix_bin_double(filename)
%#codegen
    % READ_MATRIX_BIN Reads a 2D int64 matrix from a binary file.
    % 
    % Inputs:
    %   filename - String, name of the binary file (without extension)
    %
    % Output:
    %   matrix - The reconstructed 2D int64 matrix

    binfile = [filename, '.bin']; % Binary file containing matrix data
    metafile = [filename, '.meta']; % Metadata file with matrix dimensions
    
    if ~isfile(binfile) || ~isfile(metafile)
        error('Binary or metadata file is missing.');
    end
    
    % Read the stored matrix dimensions as int64
    fid = fopen(metafile, 'rb');
    dims = fread(fid, 2, 'int64'); % Read row and column size
    fclose(fid);
    
    numRows = dims(1);
    numCols = dims(2);
    
    % Read matrix from binary file
    fid = fopen(binfile, 'rb');
    if fid == -1
        error('Cannot open binary file for reading.');
    end
    matrix = fread(fid, [numRows, numCols], 'double'); % Read as float64
    fclose(fid);
end
