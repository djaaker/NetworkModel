classdef FIFOSpikes
    properties
        buffer      % 2D matrix to store columns
        head = 1    % Index of the oldest column
        tail = 1    % Index for inserting new column
        capacity    % Maximum number of columns
        isFull = false  % Indicates if the buffer is full
    end
    methods
        function obj = FIFOSpikes(numRows, maxCols)
            obj.capacity = maxCols;
            obj.buffer = false(numRows, maxCols); % Preallocate with NaNs
        end
        
        function obj = push(obj, colVector)
            if length(colVector) ~= size(obj.buffer, 1)
                error('Column size must match the number of rows in the buffer');
            end
            
            obj.buffer(:, obj.tail) = colVector;
            
            if obj.isFull
                obj.head = mod(obj.head, obj.capacity) + 1; % Move head forward
            end
            
            obj.tail = mod(obj.tail, obj.capacity) + 1;
            obj.isFull = obj.tail == obj.head;
        end
        
        function [obj, col] = pop(obj)
            if obj.head == obj.tail && ~obj.isFull
                error('Buffer is empty');
            end
            col = obj.buffer(:, obj.head);
            obj.buffer(:, obj.head) = nan; % Optional: Clear old column
            
            obj.head = mod(obj.head, obj.capacity) + 1;
            obj.isFull = false;
        end
        
        function print(obj)
            disp('FIFO Buffer:');
            disp(obj.buffer);
        end

        function colVector = get_behind(fifo, offsets)
            % Get the buffer size
            numRows = size(fifo.buffer, 1);
            
            % Initialize the output vector
            colVector = nan(numRows, 1);
            
            % Compute circular indices
            for i = 1:numRows
                index = mod(fifo.head - 1 - offsets(i), fifo.capacity) + 1; % Circular indexing
                colVector(i) = fifo.buffer(i, index); % Extract value from computed index
            end
        end
        function obj = push_multiple(obj, colMatrix)
            [numRows, numCols] = size(colMatrix);
            
            if numRows ~= size(obj.buffer, 1)
                error('Column matrix row size must match the FIFO buffer row size.');
            end
            
            for j = 1:numCols
                obj = obj.push(colMatrix(:, j)); % Push each column one by one
            end
        end
    end
end

% % Use
% fifo = FIFO2D_Columns(4, 5); % 4 rows, max 5 columns
% 
% % Push some columns into the buffer
% fifo = fifo.push([1; 2; 3; 4]); 
% fifo = fifo.push([5; 6; 7; 8]);
% fifo = fifo.push([9; 10; 11; 12]);
% fifo = fifo.push([13; 14; 15; 16]);
% fifo = fifo.push([17; 18; 19; 20]);
% 
% fifo.display(); % Show the buffer
% 
% % Get a column where each row's value is 2 steps behind the head
% offsets = [2; 3; 6; 1]; % All rows look 2 columns back
% colVector = get_behind(fifo, offsets)
% 
% % Output should be the column from 2 steps behind in the circular buffer

