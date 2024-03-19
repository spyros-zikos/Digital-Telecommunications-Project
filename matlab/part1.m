I = imread('parrot.png');
format short;

% count each possible pixel value 0-255
I_symbols = zeros(1,256);
for i = I
    for j = i'
        I_symbols(1,double(j)+1) = I_symbols(1,double(j)+1) + 1;
    end
end

% actual pixel values and their probabilities
symbols = zeros(1,nnz(I_symbols));
probabilities = zeros(1,nnz(I_symbols));
k=1;   % 1-16
for i = linspace(1,256,256)
    if (I_symbols(1,i) ~= 0)
        symbols(1,k) = i-1;
        % probability = symbol appearance / total appearances
        probabilities(1,k) = I_symbols(1,i)/(200*150);
        k = k + 1;
    end
end

[dict, avglen] = huffmandict(symbols, probabilities);

hf = 0;
for p = probabilities
    hf = hf - (p*log2(p));
end

fprintf("\tentropy   length    performance\n")
disp([hf, avglen, hf/avglen])

% ########### 2nd order ############
% count each possible combination of 2 pixel values 0 - 65535
I_symbols2 = zeros(1,65536);
subsum = 0;
even = 0;
for i = I
    for j = i'
        if (even ~= 0)
            I_symbols2(1,double(j)+subsum+1) = ...
            I_symbols2(1,double(j)+subsum+1) + 1;
            even = 0;
        else
            subsum = double(j)*256;
            even = 1;
        end
    end
end

% actual pixel values and their probabilities
symbols2 = zeros(1,nnz(I_symbols2));
probabilities2 = zeros(1,nnz(I_symbols2));
k=1;   % 1-186
for i = linspace(1,65536,65536)
    if (I_symbols2(1,i) ~= 0)
        symbols2(1,k) = i-1;
        % probability = symbol appearance / total appearances
        probabilities2(1,k) = I_symbols2(1,i)/(100*150);
        k = k + 1;
    end
end

[dict2, avglen2] = huffmandict(symbols2, probabilities2);

hf2 = 0;
for p = probabilities2
    hf2 = hf2 - (p*log2(p));
end

fprintf("\tentropy   length    performance\n")
disp([hf2, avglen2, hf2/avglen2])

% ########### codify source ############
enc = huffmanenco(I(:)',dict);
deco = huffmandeco(enc, dict);

if (deco' == I(:))
    fprintf("\tDecoding the encoded source we get the initial source.\n")
end


bI = reshape((dec2bin(typecast(I(:),'uint8'),4)-'0').',1,[]);
J = size(enc,2)/size(bI,2);
fprintf("\tJ = %f\n\n", J)

% ########### channel ############
x=enc';

zeros_x = 0;
zeros_y_correct = 0;

many_times = 100;   % ideal: 10000+
for k=(1:many_times)
    y=binary_symmetric_channel(x);
    for i=(1:size(x))
        if x(i) == 0
            zeros_x = zeros_x + 1;
            if y(i) == 0
                zeros_y_correct = zeros_y_correct + 1;
            end
        end
    end
end
p_correct = round(zeros_y_correct/zeros_x,2);
fprintf("\tProbability p: %.2f\n", p_correct)

x_size = size(x);
p_x0 = (zeros_x/many_times)/x_size(1);
fprintf("\tProbability of x0=0: %.4f\n", p_x0)
