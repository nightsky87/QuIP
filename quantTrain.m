function quantTrain(numPatches)
    if nargin == 0
        numPatches = 1000000;
    end
    
    % Collect training sets of patch triples for training
    P = collectPatchTriples(numPatches);
    
    % Generate the forward DCT matrix for the luma channel
    [Df,~] = dctMatrix(8);
    
    % Apply the DCT to the luma patches
    TL = floor(Df' * P{1}(1:64,:));
    TR = floor(Df' * P{1}(129:end,:));

    % Store only the DC difference
    TR(1,:) = TR(1,:) - TL(1,:);

    % Store the coefficients
    P{1} = uint16(abs(TR));
    
    % Generate the forward DCT matrix for the luma channel
    [Df,~] = dctMatrix(4);
    
    % Generate a downsampling operator for the chroma channel
    Ds = zeros(16,64);
    a = 1;
    for j = 1:2:8
        for i = 1:2:8
            T = zeros(8);
            T(i:i+1,j:j+1) = 0.25;
            Ds(a,:) = T(:);
            a = a + 1;
        end
    end

    % Process the chroma channels
    for i = 2:3
        % Apply the DCT to the luma patches
        TL = floor(Df' * Ds * P{i}(1:64,:));
        TR = floor(Df' * Ds * P{i}(129:end,:));

        % Store only the DC difference
        TR(1,:) = TR(1,:) - TL(1,:);

        % Store the coefficients
        P{i} = uint16(abs(TR));
    end
    
    % Generate the rate-distortion model of the coefficients and optimize
    Q = cell(1,3);
    for i = 1:3
        % Determine the size of the patch
        numDims = size(P{i},1);
        
        % Find the rate-distortion response
        [LD,QC,TR,TD] = rdTables(P{i});
        
        % Determine the maximum rate from the table
        r = find(LD(end,:) ~= 1e99,1,'last');
        
        % Determine the optimal quantizer
        Qb = zeros(numDims,1);
        Qb(end) = QC(end,r);
        for j = numDims-1:-1:1
            r = r - TR(j+1,Qb(j+1));
            Qb(j) = QC(j,r);
        end
        
        % Determine the closest known rate to 0.5 bpp
        vr = find(LD(end,:) ~= 1e99);
        r = vr(find(vr >= 0.5*128*numDims,1,'first'));
        
        % Determine the optimal quantizer
        Qa = zeros(numDims,1);
        Qa(end) = QC(end,r);
        for j = numDims-1:-1:1
            r = r - TR(j+1,Qa(j+1));
            Qa(j) = QC(j,r);
        end
        
        % Determine the closest known rate to 0.05 bpp
        r = vr(find(vr >= 0.05*128*numDims,1,'first'));
        
        % Determine the optimal quantizer
        Qw = zeros(numDims,1);
        Qw(end) = QC(end,r);
        for j = numDims-1:-1:1
            r = r - TR(j+1,Qw(j+1));
            Qw(j) = QC(j,r);
        end
        
        % Combine the three quantization points
        Q{i} = [Qw Qa Qb];
    end
    
    % Find the average of the chroma quantization tables
    Q{2} = round((Q{2} + Q{3}) / 2);
    Q(3) = [];
    
    % Save the quantization tables to a file
    save('matrices/quantTables.mat','Q');
end

%% Sub-function for collecting training patches
function P = collectPatchTriples(numPatches)
    % Divide the number of patches between the 24 Kodak images
    N = floor(numPatches / 24) * ones(1,24);
    N(1:numPatches-sum(N)) = N(1:numPatches-sum(N)) + 1;
    
    % Process each Kodak image
    P = cell(3,24);
    for i = 1:24
        % Load the image to memory
        X = imread(sprintf('train/kodim%02d.png',i));
        
        % Map to YCbCr
        X = double(rgb2ycbcr(X));
        
        % Collect random triple-width patches from each channel
        P{1,i} = im2colrand(X(:,:,1),[8 24],N(i));
        P{2,i} = im2colrand(X(:,:,2),[8 24],N(i));
        P{3,i} = im2colrand(X(:,:,3),[8 24],N(i));
    end
    
    % Merge all sets
    P{1} = single(cell2mat(P(1,:)));
    P{2} = single(cell2mat(P(2,:)));
    P{3} = single(cell2mat(P(3,:)));
    P = P(:,1);
end
