% Change arguments according to do_everything documentation (comments)
do_everything(2,1)


function do_everything(question_num, source_mode)
% question_num: 2 for question 2, 3 for question 3 and 4 for question 4
% source_mode: 1 for the initial source, 2 for first 5000 samples, 3 for 
% second 5000 samples, 4 for third 5000 samples, 5 for last 5000 samples
    S = load('source.mat').t;
    
    if source_mode==2
        S = S(1:5000);
    elseif source_mode==3
        S = S(5001:10000);
    elseif source_mode==4
        S = S(10001:15000);
    elseif source_mode==5
        S = S(15001:20000);
    end
    s = size(S,1);

    if question_num==2
        for p=[5,10]
            figure()
            for N=[1,2,3]
                [~, ~, y] = transmitter(S, p, N);
                subplot(3,1,N)
                plot(1:s, S, '*')
                hold on
                plot(1:s, y, '*')
                hold off
                title('Initial signal and prediction errors',...
                    ['p:',int2str(p), ', N:', int2str(N)])
            end
        end
    elseif question_num==3
        pc = 1;   % point counter
        points_y = zeros(6*3, 1);
        labels = '';
        for p=5:10
            for N=1:3
                [~, ~, y] = transmitter(S, p, N);
                y = y.^2;
                y = sum(y)/size(y,1);
                
                points_y(pc)=y;
                labels = [labels 'p:' int2str(p) ', N:' int2str(N),'-'];
                pc = pc + 1;
            end
        end
        figure()
        plot(1:pc-1, points_y', "-*")
        labels = split(labels,'-');
        text(1:pc-1,points_y',labels(1:18),'VerticalAlignment', ...
            'bottom','HorizontalAlignment','right')
        title('MSE for p=5:10 and N=1:3')

        as = zeros(6,10);
        legend_txt = '';
        figure()
        for p=5:10
            [~, a, ~] = transmitter(S,p,3);
            as(p-4,1:p) = a';
            plot(1:p,a,'-x')
            legend_txt = [legend_txt,'p:',int2str(p),'-'];
            hold on
        end
        hold off
    
        legend_txt = split(legend_txt, '-');
        legend(legend_txt(1:end-1))
        title('Coefficients a, for p=5:10')
        disp('Coefficients:')
        disp(as)
    elseif question_num==4
        for p=[5,10]
            figure()            
            for N=1:3
                [yc, aq, ~] = transmitter(S, p, N);
                S_tonos = receiver(yc, aq, p);
                subplot(3,1,N)
                plot(S,'-')
                hold on
                subplot(3,1,N)
                plot(S_tonos,'-')
                hold off
                legend('Initial Signal','Reconstructed Signal')
                title(['p:' int2str(p) ', N:' int2str(N)])
            end
        end
    end
end


function [y_caps,aq,ys] = transmitter(x, p, N)
    a = coef(x,p);
    aq = zeros(size(a,1),1);
    for i=1:p
        [aq_pos,aq_centers] = my_quantizer(a(i),8,-2,2);
        aq(i) = aq_centers(aq_pos);
    end

    memory = zeros(p, 1);
    y_caps = zeros(size(x,1),1);
    ys = zeros(size(x,1),1);

    for i=1:size(x,1)

        y_tonos = prediction(memory, aq, p);

        y = x(i)-y_tonos;
        ys(i)=y;

        [y_cap_pos, y_cap_centers] = my_quantizer(y, N, -3.5, 3.5);
        y_cap = y_cap_centers(y_cap_pos);
        y_caps(i) = y_cap;

        y_cap_tonos = y_cap + y_tonos;
        
        memory = [memory(2:p); y_cap_tonos];
    end
end


function [y_cap_tonos] = receiver(y_caps, aq, p)
    y_cap_tonos = zeros(size(y_caps,1), 1);
    memory = zeros(p, 1);

    for i=1:size(y_caps,1)
        y_tonos = prediction(memory, aq, p);

        y_cap_tonos(i) = y_tonos + y_caps(i);

        memory = [memory(2:p); y_cap_tonos(i)];
    end
end


function [a,R,r] = coef(x, p)
    r = zeros(p,1);
    for i=1:p
        r(i) = autocorrelation(x,size(x,1),p,i,0);
    end

    R = zeros(p,p);
    for i=0:p-1
        for j=0:p-1
            R(i+1,j+1) = autocorrelation(x,size(x,1),p,i,j);
        end
    end
    
    a = R \ r;
end


function [Rij] = autocorrelation(x,N,p,i,j)
    Rij = 0;
    for n=p+1:N
        Rij = Rij + x(n-i)*x(n-j);
    end
    Rij = Rij/(N-p);
end


function [y_tonos] = prediction(y_cap_tonos, aq, p)
    y_tonos = 0;
    for i=1:p
        y_tonos = y_tonos + aq(i)*y_cap_tonos(p-i+1);
    end
end


function [pos, centers] = my_quantizer(yn,N,min_value,max_value)
    if yn>max_value
        yn=max_value;
    elseif yn<min_value
        yn=min_value;
    end
    
    levels = 2^N;
    delta = (max_value-min_value)/levels;
    
    centers(1) = max_value - (delta/2);
    for i=2:levels
        centers(i) = centers(i-1) - delta;
    end
    
    for i=1:size(centers')
        if (centers(i)+(delta/2)>=yn) && (centers(i)-(delta/2)<=yn)
            pos=i;
            break;
        end
    end
end
