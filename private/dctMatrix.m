%% Function for generating a DCT matrix and different scan orders
function [Df,Di] = dctMatrix(bsize)
    % Manually define the 16x16 1-D HEVC DCT matrix
    D = [64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64,   64, ...
        90,   87,   80,   70,   57,   43,   25,    9,   -9,  -25,  -43,  -57,  -70,  -80,  -87,  -90, ...
        89,   75,   50,   18,  -18,  -50,  -75,  -89,  -89,  -75,  -50,  -18,   18,   50,   75,   89, ...
        87,   57,    9,  -43,  -80,  -90,  -70,  -25,   25,   70,   90,   80,   43,   -9,  -57,  -87, ...
        83,   36,  -36,  -83,  -83,  -36,   36,   83,   83,   36,  -36,  -83,  -83,  -36,   36,   83, ...
        80,    9,  -70,  -87,  -25,   57,   90,   43,  -43,  -90,  -57,   25,   87,   70,   -9,  -80, ...
        75,  -18,  -89,  -50,   50,   89,   18,  -75,  -75,   18,   89,   50,  -50,  -89,  -18,   75, ...
        70,  -43,  -87,    9,   90,   25,  -80,  -57,   57,   80,  -25,  -90,   -9,   87,   43,  -70, ...
        64,  -64,  -64,   64,   64,  -64,  -64,   64,   64,  -64,  -64,   64,   64,  -64,  -64,   64, ...
        57,  -80,  -25,   90,   -9,  -87,   43,   70,  -70,  -43,   87,    9,  -90,   25,   80,  -57, ...
        50,  -89,   18,   75,  -75,  -18,   89,  -50,  -50,   89,  -18,  -75,   75,   18,  -89,   50, ...
        43,  -90,   57,   25,  -87,   70,    9,  -80,   80,   -9,  -70,   87,  -25,  -57,   90,  -43, ...
        36,  -83,   83,  -36,  -36,   83,  -83,   36,   36,  -83,   83,  -36,  -36,   83,  -83,   36, ...
        25,  -70,   90,  -80,   43,    9,  -57,   87,  -87,   57,   -9,  -43,   80,  -90,   70,  -25, ...
        18,  -50,   75,  -89,   89,  -75,   50,  -18,  -18,   50,  -75,   89,  -89,   75,  -50,   18, ...
        9,  -25,   43,  -57,   70,  -80,   87,  -90,   90,  -87,   80,  -70,   57,  -43,   25,   -9];
    
    % Reshape the matrix
    D = reshape(D,[16 16]);
    
    % Downsample to the correct size
    if bsize == 4
        D = D(1:4,1:4:end);
    elseif bsize == 8
        D = D(1:8,1:2:end);
    end
    
    % Generate the 2D DCT matrix
    D2 = zeros(bsize*bsize,'single');
    a = 1;
    for j = 1:bsize
        for i = 1:bsize
            T = repmat(D(:,i),[1 bsize]) .* repmat(D(:,j)',[bsize 1]);
            D2(:,a) = reshape(T,[],1);
            a = a + 1;
        end
    end
    D = D2;
    
    % Scale the dictionary elements
    if bsize == 4
        Df = D / 512;
        Di = D / 128 / 4096;
    elseif bsize == 8
        Df = D / 2048;
        Di = D / 128 / 4096;
    elseif bsize == 16
        Df = D / 8192;
        Di = D / 128 / 4096;
    end
    
    % Define the horizontal raster scan
    indV = 1:bsize*bsize;
    
    % Define the vertical raster scan
    indH = [];
    for i = 1:bsize
        indH = [indH (i:bsize:bsize*bsize)];
    end
    
    % Define the diagonal scan
    indD = [];
    for i = 1:bsize
        indD = [indD (i:bsize-1:(i-1)*bsize+1)];
    end
    for i = 2:bsize
        indD = [indD (i*bsize:bsize-1:bsize*bsize-bsize+i)];
    end
    
    % Reorder the DCT coefficients
    Df = Df(:,indD);
    Di = Di(:,indD);
end