function matrix = distribute_numbers_square(n)
    % DISTRIBUTE_NUMBERS_SQUARE Distributes consecutive numbers into a square matrix
    % Each number occupies a square patch of 10x10
    
    % Determine the size of the square matrix
    patch_size = 10; % Each number occupies a 10x10 patch
    grid_size = ceil(sqrt(n)); % Minimum grid size needed
    matrix_size = grid_size * patch_size; % Full matrix size
    
    % Initialize the matrix
    matrix = zeros(matrix_size);
    
    % Fill the matrix with consecutive numbers
    number = 1;
    for i = 1:patch_size:matrix_size
        for j = 1:patch_size:matrix_size
            if number <= n
                matrix(i:i+patch_size-1, j:j+patch_size-1) = number;
                number = number + 1;
            end
        end
    end
    
    % Trim excess parts if the last number is reached before filling the grid
    matrix = matrix(1:grid_size*patch_size, 1:grid_size*patch_size);
end
