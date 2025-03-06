function save_matrix_bin(filename, matrix, mode)
%#codegen
    % SAVE_MATRIX_BIN Saves or appends a 2D int64 matrix to a binary file.
    % 
    % Inputs:
    %   filename - String, name of the binary file (without extension)
    %   matrix   - 2D int64 matrix to be saved/appended
    %   mode     - 'w' to write a new file, 'a' to append columns
    
    binfile = [filename, '.bin']; % Binary file to store matrix data
    metafile = [filename, '.meta']; % Metadata file to store matrix dimensions
    
    if nargin < 3
        mode = 'w'; % Default to write mode
    end
    
    if mode == 'w'
        % Overwrite and save the matrix
        fid = fopen(binfile, 'wb');
        if fid == -1
            error('Cannot open binary file for writing.');
        end
        fwrite(fid, matrix, 'double'); % Store as int64
        fclose(fid);
        
        % Save dimensions separately as int64
        fid = fopen(metafile, 'wb');
        if fid == -1
            error('Cannot open metadata file for writing.');
        end
        fwrite(fid, int64([size(matrix, 1), size(matrix, 2)]), 'int64');
        fclose(fid);
        
    elseif mode == 'a'
        % Ensure the metadata file exists
        if ~isfile(metafile)
            error('Metadata file missing. Use mode "w" to create a new file first.');
        end
        
        % Read stored matrix dimensions as int64
        fid = fopen(metafile, 'rb');
        dims = fread(fid, 2, 'int64');
        fclose(fid);
        
        numRows = dims(1);
        numCols = dims(2);
        
        % Check row consistency
        if size(matrix, 1) ~= numRows
            error('New data must have the same number of rows as the existing matrix.');
        end
        
        % Append new columns to the binary file
        fid = fopen(binfile, 'ab');
        if fid == -1
            error('Cannot open binary file for appending.');
        end
        fwrite(fid, matrix, 'double'); % Append in int64 format
        fclose(fid);
        
        % Update the metadata file with the new column count
        fid = fopen(metafile, 'wb');
        if fid == -1
            error('Cannot open metadata file for updating.');
        end
        fwrite(fid, int64([numRows, numCols + size(matrix, 2)]), 'int64');
        fclose(fid);
    else
        error('Invalid mode. Use "w" for write or "a" for append.');
    end
end
