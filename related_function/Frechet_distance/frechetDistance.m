function [frechetDist, couplingSequence, frechetMatrix, distanceMatrix] = frechetDistance(curve1, curve2, distanceMatrixFunction, endpointFlexibilityFlag)
% Calculates the discrete Frechet distance between curve1 and curve2.  
% 
% Simple usage:
% frechetDist = frechetDistance(curve1, curve2)
% 
% Returned in frechetDist is the distcrete Frechet distance between the
% curves supplied in curve1 and curve2. The curves should have one point
% per row, and may be of arbitrary (but matching) dimensionality. In simple
% usage, the distance is the n-dimensional Euclidean distance, and endpoint
% flexibility is not allowed (i.e. the sequence of point couplings must
% include pairing the starting points and pairing the ending points of each
% curve). 
% 
% Complex usage:
% [frechetDist, couplingSequence, frechetMatrix, distanceMatrix] = ...
%    frechetDistance(curve1, curve2, distanceMatrixFunction, endpointFlexibilityFlag)
% 
% The output frechetDist and the inputs curve1 and curve2 are the same as
% in simple usage.  
% 
% Optional input distanceMatrixFunction should be a function handle for a
% function which takes in curve1 and curve2 as inputs and which returns the
% pairwise distance matrix between the points of curve1 and the points of
% curve2. That is, if dm = distanceMatrixFunction(curve1, curve2), then
% dm(i,j) is the distance between the point curve1(i,:) and the point
% curve2(j,:).  This allows the user to specify any arbitrary desired
% distance metric calculation, and would even allow a distance to be
% specified between points of differing dimensionalities (if such a thing
% were desired). If this input is empty ([]) or the number of arguments is
% less than 3, distanceMatrixFunction defaults to n-dimensional Euclidean
% distance (i.e. the L2 norm between points).
% 
% Optional input endpointFlexibilityFlag should be either true or false. If
% false, empty, or there are fewer than 4 input arguments, endpoint
% flexibility is not allowed, meaning that the sequence of point couplings
% must include pairing the starting points and ending points of each curve
% (which is standard for typical Frechet distance calculations). However,
% for some applications, this behavior is misleading about the
% dissimilarity between curves. For example, consider two periodic curves
% defined as follows:
% x1 = linspace(0,2*pi,500)';
% y1 = sin(10.*x1);
% curve1 = [x1, y1];
% x2 = x1;
% y2 = sin(10.*(x2-0.1));
% curve2 = [x2, y2];
% Intuitively, the Frechet distance between these curves should be 0.1,
% since curve2 is just curve1 shifted to the right by 0.1.  However, the
% Frechet distance between these two curves comes out to 0.8415 rather than
% 0.1, because the standard processing forces the starting points to be
% paired together, and these points are 0.8415 apart.  If the curves
% continued infinitely, then every point could be paired with one 0.1 away,
% but since the curves are out of phase and we are looking at the same
% range of x values, the endpoints are typically more than 0.1 apart. The
% endpointFlexibilityFlag allows this type of issue to be mitigated by
% relaxing the requirement that both starting points (and ending points)
% must be paired together.  Instead, the distance between the first point
% of curve1 and the closest point on curve2 to it is found, and the
% distance between the first point of curve2 and the closest point on
% curve1 is found, and whichever of these is shorter forms the first link
% in the coupling sequence. Essentially, this allows ignoring the early
% part of one (but not both) of the two curves when calculating the Frechet
% distance.  The same process is used at the other end of the curves,
% allowing the tail end of one of the curves to be ignored if there is a
% closer link between one of the ends and an earlier point on the other
% curve. If we run this function again on the curves defined above and
% allow endpoint flexibility, we find that the calculated Frechet distance
% falls to essentially 0.1 as expected (it's slightly greater because of
% the discrete sampling of the curves). It is recommended to think
% carefully about your data and whether it would be appropriate to allow
% endpoint flexibility or not. The easiest way to see the impact is to use
% plotCouplingSeq.m with a coupling sequence generated without endpoint
% flexibility and then with a coupling sequence generated with endpoint
% flexibility and compare them.
% 
% Output couplingSequence is an nx2 matrix of point indices representing
% one possible sequence of curve1 to curve2 point couplings with the
% calculated Frechet distence;  that is, the distance between the indicated
% pairs of points is equal to the Frechet distance for at least one
% pairing, and less than or equal to the Frechet distance for all other
% pairings. The first column of couplingSequence is a row index into
% curve1, and the second column of couplingSequence is a row index into
% curve2. Each column vector has the property that the next entry is either
% equal to or one greater than the previous entry, guaranteeing that curve
% points are visited in sequential order. If endpointFlexibility is not
% allowed, the first row of couplingSequence will always be [1,1] and the
% last row will be [n1,n2] if there are n1 points in curve1 and n2 points
% in curve2. If endpointFlexibility is allowed, then the first row will
% include at least one 1 and the last row will include at least one of n1
% or n2. Note that there may be many possible coupling sequences with the
% same Frechet distance, since there may be many options for couplings
% which are not minimal but do not exceed the Frechet distance. This
% function tries to build a "good" coupling sequence by preferring
% couplings which are shorter distance to those which are longer, but even
% these are not guaranteed to be unique. The returned couplingSequence
% should be considered an illustrative example of a possible coupling which
% satisfies the calculated frechetDist. 
% 
% Output frechetMatrix is, if endpointFlexibility is not allowed, an n1 by
% n2 matrix where the frechetMatrix(i,j) is the Frechet distance between
% the curves curve1(1:i,:) and curve2(1:j,:). Thus, the Frechet distance
% between the full curve1 and curve2 is frechetMatrix(n1,n2). The
% frechetMatrix is iteratively or recursively constructed starting from the
% pairwise distance matrix by considering the minimal longest link to get
% to a certain pairing of points. If endpointFlexibility is allowed, then
% frechetMatrix may be smaller than n1 by n2, because ignored points near
% either end will not be included.  For example, if the first 3 points of
% curve1 are ignored and the last 2 points of curve2 are ignored, then the
% output frechetMatrix will be (n1-3) by (n2-2). The calculated
% frechetDist will still be the last-row, last-column entry of the
% frechetMatrix. 
% 
% Output distanceMatrix is the full pairwise distance matrix between the
% points of curve1 and curve2.  It will always be n1 by n2 regardless of
% whether endpoint flexibility is allowed. 



%% Construct Pairwise Distance Matrix
if nargin<3 || isempty(distanceMatrixFunction)
    % Default to normal euclidean distance between points in n dimensions
    distanceMatrixFunction = @defaultDistanceMatrixFunction;
end

distanceMatrix = distanceMatrixFunction(curve1, curve2);

%% Handle endpoint flexibility
if nargin<4 || isempty(endpointFlexibilityFlag)
    endpointFlexibilityFlag = false;
end
if endpointFlexibilityFlag
    % Allow one of the curves to start later than the first point if that
    % allows a better match
    [minVal1,minIdx1] = min(distanceMatrix(:,1));
    [minVal2,minIdx2] = min(distanceMatrix(1,:));
    if minVal1 <= minVal2
        startPoint1 = minIdx1;
        startPoint2 = 1;
    else
        startPoint1 = 1;
        startPoint2 = minIdx2;
    end
    % Do the same for the ends of the curves
    [minEndVal1, minEndIdx1] = min(distanceMatrix(:,end));
    [minEndVal2, minEndIdx2] = min(distanceMatrix(end,:));
    if minEndVal1 <= minEndVal2
        endPoint1 = minEndIdx1;
        endPoint2 = size(distanceMatrix,2);
    else
        endPoint1 = size(distanceMatrix,1);
        endPoint2 = minEndIdx2;
    end
else
    % Force use of full curves
    startPoint1 = 1;
    startPoint2 = 1;
    endPoint1 = size(distanceMatrix,1);
    endPoint2 = size(distanceMatrix,2);
end
trimmedDistanceMatrix = distanceMatrix(startPoint1:endPoint1, startPoint2:endPoint2);

%% Evolve the Frechet Matrix
% In the end, this matrix has in each element the Frechet distance
% necessary to reach that pair of points.
frechetMatrix = evolveFrechetMatrix3(trimmedDistanceMatrix);


%% Get the Frechet Distance of the full pair of curves
frechetDist = frechetMatrix(end,end);

%% Build a coupling sequence
couplingSequence = buildCouplingSequence(frechetMatrix, trimmedDistanceMatrix);
% If there were any trimmed starting points, add them back to the index
% values (couplingSequence values are indices into the trimmed distance
% matrix, for these to be indices into the original curve points, we need
% to add the starting point index
% minus 1).
couplingSequence(:,1) = couplingSequence(:,1) + startPoint1 - 1;
couplingSequence(:,2) = couplingSequence(:,2) + startPoint2 - 1;



%% Subfunctions
%% Full matrix version of frechetMatrix calculation (vectorized but very slow)
    function frechetMatrix = evolveFrechetMatrix1(distanceMatrix)
        % In the end, this matrix has in each element the Frechet distance
        % necessary to reach that pair of points. This is the slowest
        % approach, taking 20 sec on a 1000x1000 distance matrix. In each
        % iteration, an 'alternative' matrix is generated which has the
        % minimum value of the left, left-up, and up neighboring values,
        % and then the frechetMatrix is updated to be the element-wise
        % maximum of either the alternative matirx or the current value.
        % This forward-propagates values to the right and down. In the
        % worst case scenario, where the maximum distance is between the
        % first two points, it will take max([n1,n2]) iterations for that
        % value to propagate to the last row and column if [n1,n2] =
        % size(distanceMatrix). However, it may take fewer iterations for
        % the matrix to reach it's final state, and if the matrix is ever
        % unchanged from one iteration to the next, then it has finished
        % evolving and computation can cease. 

        % Initialize
        frechetMatrix = distanceMatrix;
        sz = size(distanceMatrix);
        % Iterate at most the maximum dimension (the moving front of final
        % answers must advance at least one row and one column each step)
        for iter = 1:max(sz)
            frechetMatrixPrev = frechetMatrix;
            % Construct an array of alternative values (the minimum of the
            % left, up, and left-up neighbors)
            alt = zeros(sz);
            leftEntries = frechetMatrix(:,1:end-1);
            upEntries = frechetMatrix(1:end-1,:);
            leftUpEntries = frechetMatrix(1:end-1,1:end-1);
            % Handle the bulk of the alternative matrix
            alt(2:end,2:end) = min(cat(3, leftEntries(2:end,:), upEntries(:,2:end), leftUpEntries), [], 3);
            % Handle the upper and left edges
            alt(1,2:end) = leftEntries(1,:);
            alt(2:end,1) = upEntries(:,1);
            alt(1,1) = distanceMatrix(1,1);
            % Update frechetMatrix with the maximum of the alt matrix or the
            % distance matrix
            frechetMatrix = max(distanceMatrix, alt);
            % Check if it changed anything
            if all(frechetMatrix(:) == frechetMatrixPrev (:))
                % No longer evolving, we can stop calculating
                %disp(sprintf('Breaking at step %i', iter))
                break
            end
        end
    end
%% Forward neighbor version of frechetMatrix calculation (fast)
    function frechetMatrix = evolveFrechetMatrix2(distanceMatrix)
        % In this approach, instead of calculating the full alternative
        % matrix in each round, only the neighbors of values which were
        % changed in the previous round are checked.  Updates are continued
        % until there is a round with no changes. This is much faster (>10x
        % faster) than the full matrix version, taking 1.75 seconds on a
        % sample 1000x1000 distance matrix.

        frechetMatrix = distanceMatrix; % initialize

        %%% Edge rows
        % The top row has only left neighbors and the left column only has
        % top neighbors, so we can quickly handle those entries and evolve
        % them to their final values.  Then, all remaining entries have
        % three relevant neighbors and we won't have to handle edge cases
        % anymore.
        frechetMatrixTopRow = frechetMatrix(1,:);
        for colIdx = 2:length(frechetMatrixTopRow)
            frechetMatrixTopRow(colIdx) = max(frechetMatrixTopRow(colIdx), frechetMatrixTopRow(colIdx-1));
        end
        frechetMatrixLeftCol = frechetMatrix(:,1);
        for rowIdx = 2:length(frechetMatrixLeftCol)
            frechetMatrixLeftCol(rowIdx) = max(frechetMatrixLeftCol(rowIdx), frechetMatrixLeftCol(rowIdx-1));
        end

        frechetMatrix(1,:) = frechetMatrixTopRow;
        frechetMatrix(:,1) = frechetMatrixLeftCol;

        %%% First update
        % For the first round, we need to check everyone's neighbors, so
        % for this case we might as well do it in matrix form
        updateMask = false(size(frechetMatrix));
        leftEntries = frechetMatrix(2:end,1:end-1);
        upEntries = frechetMatrix(1:end-1,2:end);
        leftUpEntries = frechetMatrix(1:end-1,1:end-1);
        minNeightborValue = min(cat(3, leftEntries, upEntries, leftUpEntries),[],3);
        updateMask(2:end,2:end) = frechetMatrix(2:end,2:end) < minNeightborValue;
        frechetMatrix(updateMask) = minNeightborValue(updateMask(2:end,2:end));

        %%% Iterative updates
        % After the first round, we only need to check on neighbors of changed
        % entries.  This should be faster than recalculating everything
        edgeMask = false(size(updateMask));
        edgeMask(1,:) = true;
        edgeMask(:,1) = true;

        leftOffset = -1*size(frechetMatrix,1);
        upOffset = -1;
        leftUpOffset = leftOffset+upOffset;
        while sum(updateMask(:))>0
            % Update candidates are those to the right, below, or below-right of
            % updated values (and which are not in the top row or left column)
            updateCandidatesMask = ( ...
                circshift(updateMask,1,2) | ... % right
                circshift(updateMask,1,1) | ... % below
                circshift(circshift(updateMask,1,1),1,2) ... % below-right
                ) ...
                & ~edgeMask; % and not on the upper or left edge
            uc = find(updateCandidatesMask);

            currentVals = frechetMatrix(uc);
            leftVals = frechetMatrix(uc+leftOffset);
            upVals = frechetMatrix(uc+upOffset);
            leftUpVals = frechetMatrix(uc+leftUpOffset);
            altVals = min([leftVals, upVals, leftUpVals],[],2);
            toUpdate = uc(altVals > currentVals);
            newUpdateMask = false(size(frechetMatrix));
            newUpdateMask(toUpdate) = true;
            % Carry out the update
            frechetMatrix(toUpdate) = altVals(altVals>currentVals);
            updateMask = newUpdateMask;
        end % end of while
        % frechetMatrix is complete after exiting the while loop
    end
%% Recursive version of frechetMatrix calculation (fastest)
    function frechetMatrix = evolveFrechetMatrix3(distanceMatrix)
        % Recursive version.  This is the fastest approach, but one I found
        % difficult to understand and verify at first (hence the other
        % approaches I tried).  It takes 0.86 seconds on a sample 1000x1000
        % distanceMatrix (about 1/2 the time of approach 2 and about 1/23rd
        % the time of approach 1).
        frechetMatrix = -1.*ones(size(distanceMatrix));
        function fmij = c(i,j)
            if frechetMatrix(i,j) > -1
                % If entry is already filled in, just look it up
                fmij = frechetMatrix(i,j);
            elseif i==1 && j==1
                % Entry 1,1 is always just the distance between first
                % points
                fmij = distanceMatrix(1,1);
            elseif i>1 && j==1
                % Left column entries are the larger of the distance matrix
                % entry or the upper neighbor
                fmij = max( c(i-1,1), distanceMatrix(i,1) );
            elseif i==1 && j>1
                % Top row entries are the larger of the distance matrix
                % entry or the left neighbor
                fmij = max( c(i,j-1), distanceMatrix(1,j) );
            elseif i>1 && j>1
                % All other entries are the larger of the distance matrix
                % entry or (the minimum of the left, up-left, or up
                % neighbors)
                fmij = max( min( [c(i-1,j), c(i-1,j-1), c(i,j-1)]), distanceMatrix(i,j));
            else
                % Shouldn't be able to reach here with valid i and j values
                error('i and j must be valid indices into the pairwise distance matrix!')
            end
            frechetMatrix(i,j) = fmij;
        end
        % Trigger recursive filling in to the final entry
        c(size(frechetMatrix,1), size(frechetMatrix,2));
    end
%% Default distance matrix calculation function (Euclidean)
    function distanceMatrix = defaultDistanceMatrixFunction(curve1, curve2)
        % This function calculates the pairwise distance matrix function
        % between all points on two curves, i.e., distanceMatrix(i,j) is the
        % distance between the point curve1(i,:) and the point curve2(j,:).
        % The reported distance is the Euclidean distance between points in
        % n-dimensions.
        [numPoints1, numDims1] = size(curve1);
        [numPoints2, numDims2] = size(curve2);
        if numDims1 ~= numDims2
            error('Can''t calculate euclidean distance between points of different dimensionality!');
        end
        curve1PointArray = repmat(permute(curve1, [1,3,2]), [1, numPoints2, 1]);
        curve2PointArray = repmat(permute(curve2, [3,1,2]), [numPoints1, 1, 1]);
        distanceMatrix = vecnorm(curve1PointArray-curve2PointArray, 2, 3);
    end
%% Build coupling sequence
    function couplingSequence = buildCouplingSequence(frechetMatrix, trimmedDistanceMatrix)
        % This is a sequence of point pairings of which the maximum
        % distance between any paired point is the Frechet distance.  These
        % are typically not unique because there are often many such paths
        % with the same maximum length. The approach here tries to build a
        % "good" path which has the minimizes the maximum pair distance and
        % then, within those options, also minimizes the pair distance at
        % each coupling step, and then, within those options, prefers paths
        % through the pairing matrix which proceed most directly to the
        % destination corner.  This path is also not guaranteed to be
        % unique, but it should generally select a path which is 'nicer'
        % than alternatives, and will behave more like a dog on an elastic
        % leash where the paired points are preferred to stay closer
        % together when they can (other approaches do not care how far the
        % dog wanders until it would exceed the Frechet distance). Inputs
        % should be only the Frechet Matrix and Pairwise Distance Matrix
        % (trimmed version if endpoints are flexible), neither of which is
        % modified.
        leftIdx = 1;
        upIdx = 2;
        diagIdx = 3;
        options = [leftIdx, upIdx, diagIdx];
        [nPoints1, nPoints2] = size(trimmedDistanceMatrix);

        % Build sample coupling sequence (typically non-unique)
        couplingSequence = zeros(nPoints1 + nPoints2, 2);
        r = nPoints1;
        c = nPoints2;
        couplingSequence(1,:) = [nPoints1, nPoints2];
        stepCount = 2; % first step already entered
        while r>1 && c>1
            % Consider the three neighbors
            leftNeighbor = {r,c-1};
            upNeighbor = {r-1,c};
            leftUpNeighbor = {r-1,c-1};
            % First criterion, which has the minimum frechetDist value?
            frechetVals = [frechetMatrix(leftNeighbor{:}),...
                frechetMatrix(upNeighbor{:}), ...
                frechetMatrix(leftUpNeighbor{:})];
            minFrechet = min(frechetVals);
            minFrechetMask = frechetVals == minFrechet;
            if sum(minFrechetMask) > 1
                % If more than one have the same, then prefer the one which has a
                % smaller distMatrix value.
                distMatrixVals = [trimmedDistanceMatrix(leftNeighbor{:}), trimmedDistanceMatrix(upNeighbor{:}), trimmedDistanceMatrix(leftUpNeighbor{:})];
                minDistVal = min(distMatrixVals(minFrechetMask));
                minDistValMask = distMatrixVals(minFrechetMask)==minDistVal;
                if sum(minDistValMask) > 1
                    % If more than one have the same, prefer the one which aims more
                    % directly towards the 1,1 corner
                    if r > c+1
                        % prefer up
                        preferenceOrder = [upIdx, diagIdx, leftIdx];
                    elseif c > r+1
                        % prefer left
                        preferenceOrder = [leftIdx, diagIdx, upIdx];
                    else
                        % prefer diagonal
                        if r==c
                            % prefer up as second choice (arbitrary choice)
                            preferenceOrder = [diagIdx, upIdx, leftIdx];
                        elseif r>c
                            % prefer up as second choice
                            preferenceOrder = [diagIdx, upIdx, leftIdx];
                        else % r<c
                            % prefer left as second choice
                            preferenceOrder = [diagIdx, leftIdx, upIdx];
                        end
                    end
                    % Choose the highest preference of those with the same
                    % minimum distance value
                    frechetOptions = options(minFrechetMask);
                    distOptions = frechetOptions(minDistValMask);
                    if intersect(preferenceOrder(1), distOptions)
                        chosen = preferenceOrder(1);
                    elseif intersect(preferenceOrder(2), distOptions)
                        chosen = preferenceOrder(2);
                    else
                        error('Something went wrong. There should be at least two options at this point, and so it shouldn''t be possible to be forced to the 3rd choice!')
                    end
                else
                    % Only one choice after distance comparison
                    frechetOptions = options(minFrechetMask);
                    chosen = frechetOptions(minDistValMask);
                end
            else
                % only one choice after frechet comparison
                chosen = options(minFrechetMask);
            end
            % Update r and c and add chosen pairing to the coupling sequence
            if chosen==leftIdx
                c = c-1;
            elseif chosen==upIdx
                r = r-1;
            elseif chosen==diagIdx
                r = r-1;
                c = c-1;
            else
                error('Chosen should be 1,2, or 3, but is not.')
            end
            couplingSequence(stepCount,:) = [r,c];
            % Increment step count
            stepCount = stepCount+1;
        end % end of while
        % Trim any remaining zeros off the end of the coupling sequence
        zeroMask = couplingSequence(:,1)==0;
        couplingSequence(zeroMask,:) = [];

        % Handle any remainder of the sequence  (once we have hit the top row
        % or the left column there are no more choices, it simply has to run
        % along that row or column back to 1,1).
        if r==1 && c==1
            % The last step was a diagonal move to 1,1, we're actually done
            % already
        elseif r==1
            % Now we just have to run backwards along the top row until we get
            % back to 1,1
            cs = (c:-1:1)';
            rs = ones(size(cs));
            remainingCouplingSteps = [rs, cs];
            couplingSequence = [couplingSequence; remainingCouplingSteps];
        elseif c==1
            rs = (r:-1:1)';
            cs = ones(size(rs));
            remainingCouplingSteps = [rs, cs];
            couplingSequence = [couplingSequence; remainingCouplingSteps];
        end
        % Reverse the coupling sequence order (we built it backwards)
        couplingSequence = flipud(couplingSequence);
    end

%% Test frechetMatrix evolution functions
    function [t1,t2,t3] = testFrechetMatrixEvolutionFunctions()
        % This function was written as part of the evaluation of the three
        % methods of generating the Frechet matrix.  I had some trouble
        % understanding how the recursive method (3) was really working in a
        % way that I could verify that it was correct, so I explored other
        % methods which were clearer to me and easier to verify.  However,
        % the recursive method was the most computationally efficient, and
        % they all produce the same results, so that's the one used in the
        % main code above. 
        distMatrix = rand(1000);
        % Record timings for each method
        tic;
        frechetMatrix1 = evolveFrechetMatrix1(distMatrix);
        t1 = toc;
        tic;
        frechetMatrix2 = evolveFrechetMatrix2(distMatrix);
        t2 = toc;
        tic;
        frechetMatrix3 = evolveFrechetMatrix3(distMatrix);
        t3 = toc;
        % Ensure they all produced the same result!
        assert(all(frechetMatrix1(:)==frechetMatrix2(:)) & all(frechetMatrix1(:)==frechetMatrix3(:)));
    end
end