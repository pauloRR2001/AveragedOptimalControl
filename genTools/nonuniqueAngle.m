function val = nonuniqueAngle(array)
    
    % Number of elements in the array
    n = numel(array);
    
    % Initialize variables to track the smallest difference and its index
    minDiff = inf;
    repeatedIndex = -1;

    % Loop through each element and compare with the rest
    for i = 1:n
        for j = i+1:n
            diff = abs(array(i) - array(j));
            if diff < minDiff
                minDiff = diff;
                repeatedIndex = i; % Store the index of the repeated value
            end
        end
    end

    % Return the repeated value
    val = array(repeatedIndex);
end