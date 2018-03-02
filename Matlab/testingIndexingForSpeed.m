

myArraySize = 2.^(1 : 24);
majorRowTraversalTime = zeros(1, length(myArraySize));
majorColumnTraversalTime = zeros(1, length(myArraySize));
randomIndexingSequenceTime = zeros(1, length(myArraySize));

for nArrayElements = 1 : length(myArraySize)
    % Declare array
    myArray = zeros(ceil(sqrt(myArraySize(nArrayElements))),ceil(sqrt(myArraySize(nArrayElements))));
    
    % test row traversal speed & record
    tic();
    for iRow = 1 : size(myArray, 1)
        for jColumn = 1 : size(myArray, 2)
            myArray(iRow, jColumn) = jColumn;
        end
    end
    majorRowTraversalTime(nArrayElements) = toc();
    
    % test column traversal speed & record 
    tic();
    for iColumn = 1 : size(myArray, 1)
        for jRow = 1 : size(myArray, 2)
            myArray(iColumn, jRow) = jRow;
        end
    end
    majorColumnTraversalTime(nArrayElements) = toc();
    
    % test random indexing
    myRandomIndexingSequence =  1 + ceil((rand(1, myArraySize(nArrayElements)) .* (nArrayElements -1)));
    tic();
    for iRandomSequence = 1 : length(myRandomIndexingSequence)
        myArray(myRandomIndexingSequence) = myRandomIndexingSequence(iRandomSequence);
    end
    randomIndexingSequenceTime(nArrayElements) = toc();
end

plot(myArraySize, [majorRowTraversalTime majorColumnTraversalTime randomIndexingTraversalTime])
