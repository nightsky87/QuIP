function quantTest()
    % Load the Lena image to memory
    X = imread('test/lena_std.tif');
    
    % Map to 4:2:0 YCbCr space
    T = rgb2ycbcr(X);
    Y = single(T(:,:,1));
    Cb = single(imresize(T(:,:,2),0.5,'bicubic'));
    Cr = single(imresize(T(:,:,3),0.5,'bicubic'));
    
    % Generate DCT dictionaries
    [Df4,Di4] = dctMatrix(4);
    [Df8,Di8] = dctMatrix(8);
    
    % Vectorize the data
    DY = floor(Df8' * im2col(Y,[8 8],'distinct'));
    DCb = floor(Df4' * im2col(Cb,[4 4],'distinct'));
    DCr = floor(Df4' * im2col(Cr,[4 4],'distinct'));
    
    % Test at different quality factors
    rate = zeros(1,99);
    score = zeros(1,99);
    rate2 = zeros(1,99);
    score2 = zeros(1,99);
    for i = 1:99
        % Quantize at the given qualilty factor
        [TY,TCb,TCr] = quantCoeff(DY,DCb,DCr,i);
        [TY2,TCb2,TCr2] = quantCoeffJPEG(DY,DCb,DCr,i);
        
        % Calculate the entropy of the coefficients
        p = hist([TY(:); TCb(:); TCr(:)],0:32767);
        p = p(p > 0);
        p = p / sum(p);
        p2 = hist([TY2(:); TCb2(:); TCr2(:)],0:32767);
        p2 = p2(p2 > 0);
        p2 = p2 / sum(p2);
        
        % Reconstruct the image
        TY = col2im(floor(Di8 * TY),[8 8],[512 512],'distinct');
        TCb = col2im(floor(Di4 * TCb),[4 4],[256 256],'distinct');
        TCr = col2im(floor(Di4 * TCr),[4 4],[256 256],'distinct');
        TY2 = col2im(floor(Di8 * TY2),[8 8],[512 512],'distinct');
        TCb2 = col2im(floor(Di4 * TCb2),[4 4],[256 256],'distinct');
        TCr2 = col2im(floor(Di4 * TCr2),[4 4],[256 256],'distinct');
        
        % Combine to form the image
        T = cat(3,TY,imresize(TCb,2,'bicubic'),imresize(TCr,2,'bicubic'));
        T = ycbcr2rgb(uint8(T));
        T2 = cat(3,TY2,imresize(TCb2,2,'bicubic'),imresize(TCr2,2,'bicubic'));
        T2 = ycbcr2rgb(uint8(T2));
        
        rate(i) = -sum(p .* log2(p));
        score(i) = psnr(T,X);
        rate2(i) = -sum(p2 .* log2(p2));
        score2(i) = psnr(T2,X);
    end
    
    close all;
    plot(rate,score,'k-',rate2,score2,'k--');
end

function [Y,Cb,Cr] = quantCoeff(Y,Cb,Cr,quality)
    % Load the quantization tables
    load('matrices/quantTables.mat');
    
    % Divide the quality into two intervals
    if quality < 50
        % Calculate the weighting factor
        w = quality / 50;
        
        % Generate the quantization table
        Ql = floor(Q{1}(:,1:2) * [1-w w]' + 0.5);
        Qc = floor(Q{2}(:,1:2) * [1-w w]' + 0.5);
    else
        % Calculate the weighting factor
        w = (quality - 50) / 50;
        
        % Generate the quantization table
        Ql = floor(Q{1}(:,2:3) * [1-w w]' + 0.5);
        Qc = floor(Q{2}(:,2:3) * [1-w w]' + 0.5);
    end
    
    % Explicitly quantize the first block
    Y(:,1) = floor((Y(:,1) + floor(Ql / 2)) ./ Ql) .* Ql;
    Cb(:,1) = floor((Cb(:,1) + floor(Qc / 2)) ./ Qc) .* Qc;
    Cr(:,1) = floor((Cr(:,1) + floor(Qc / 2)) ./ Qc) .* Qc;
    
    % Initialize the DPCM for DC values
    DY = Y(1);
    DCb = Cb(1);
    DCr = Cr(1);
    
    % Sequentially quantize the blocks
    for i = 2:size(Y,2)
        % Calculate the DC difference
        TY = Y(1,i) - DY;
        TCb = Cb(1,i) - DCb;
        TCr = Cr(1,i) - DCr;
        
        % Quantize the remaining coefficients
        Y(:,i) = floor((Y(:,i) + floor(Ql / 2)) ./ Ql) .* Ql;
        Cb(:,i) = floor((Cb(:,i) + floor(Qc / 2)) ./ Qc) .* Qc;
        Cr(:,i) = floor((Cr(:,i) + floor(Qc / 2)) ./ Qc) .* Qc;
        
        % Overwrite the DC levels with the quantized DPCM value
        Y(1,i) = floor((TY + floor(Ql(1) / 2)) / Ql(1)) * Ql(1);
        Cb(1,i) = floor((TCb + floor(Qc(1) / 2)) / Qc(1)) * Qc(1);
        Cr(1,i) = floor((TCr + floor(Qc(1) / 2)) / Qc(1)) * Qc(1);
        
        % Update the DPCM values
        DY = DY + Y(1,i);
        DCb = DCb + Cb(1,i);
        DCr = DCr + Cr(1,i);
    end
    
    % Restore the DC values
    Y(1,:) = cumsum(Y(1,:));
    Cb(1,:) = cumsum(Cb(1,:));
    Cr(1,:) = cumsum(Cr(1,:));
end

function [Y,Cb,Cr] = quantCoeffJPEG(Y,Cb,Cr,quality)
    % Define the JPEG quantization table
    Ql = [16 12 11 14 12 10 14 13 14 16 18 17 16 19 24 ...
          24 22 22 24 26 40 49 35 37 29 40 58 51 ...
          72 64 55 56 51 57 60 61 92 78 64 68 87 69 55 ...
          95 87 81 109 80 56 98 103 104 103 62 ...
          112 121 113 77 100 120 92 103 101 99]';
      
    Qc = [17 18 18 24 21 24 47 26 26 47 66 56 66 99 99 99]';
    
    % Divide the quality into two intervals
    if quality < 50
        % Calculate the weighting factor
        w = 5000 / quality;
    else
        % Calculate the weighting factor
        w = 200 - 2 * quality;
    end
    
    % Scale the quantization table
    Ql = floor((128 * w * Ql + 50) / 100);
    Qc = floor((128 * w * Qc + 50) / 100);
    
    % Explicitly quantize the first block
    Y(:,1) = floor((Y(:,1) + floor(Ql / 2)) ./ Ql) .* Ql;
    Cb(:,1) = floor((Cb(:,1) + floor(Qc / 2)) ./ Qc) .* Qc;
    Cr(:,1) = floor((Cr(:,1) + floor(Qc / 2)) ./ Qc) .* Qc;
    
    % Initialize the DPCM for DC values
    DY = Y(1);
    DCb = Cb(1);
    DCr = Cr(1);
    
    % Sequentially quantize the blocks
    for i = 2:size(Y,2)
        % Calculate the DC difference
        TY = Y(1,i) - DY;
        TCb = Cb(1,i) - DCb;
        TCr = Cr(1,i) - DCr;
        
        % Quantize the remaining coefficients
        Y(:,i) = floor((Y(:,i) + floor(Ql / 2)) ./ Ql) .* Ql;
        Cb(:,i) = floor((Cb(:,i) + floor(Qc / 2)) ./ Qc) .* Qc;
        Cr(:,i) = floor((Cr(:,i) + floor(Qc / 2)) ./ Qc) .* Qc;
        
        % Overwrite the DC levels with the quantized DPCM value
        Y(1,i) = floor((TY + floor(Ql(1) / 2)) / Ql(1)) * Ql(1);
        Cb(1,i) = floor((TCb + floor(Qc(1) / 2)) / Qc(1)) * Qc(1);
        Cr(1,i) = floor((TCr + floor(Qc(1) / 2)) / Qc(1)) * Qc(1);
        
        % Update the DPCM values
        DY = DY + Y(1,i);
        DCb = DCb + Cb(1,i);
        DCr = DCr + Cr(1,i);
    end
    
    % Restore the DC values
    Y(1,:) = cumsum(Y(1,:));
    Cb(1,:) = cumsum(Cb(1,:));
    Cr(1,:) = cumsum(Cr(1,:));
end
