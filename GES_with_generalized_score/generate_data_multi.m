clear all, clc,close all
% load graph_structures2 % : G_record
rng(10)
N = 5; % N nodes
G = zeros(N, N);
G(1,2)=1; G(1,4)=1; G(2,3)=1; G(2,4)=1; G(3,4)=1; G(4,5)=1;
T = 1000;
count = 0;
func = [1,2,3]; % 1: linear, 2: sinc, 3: cos, 4: tanh
noise = [1,2]; % 1: Gaussian, 2: uniform
dmax = 2; % maximum dimension of a variable
for trial = 1:500
    X = [];
    for j = 1:N
        if(j==1)
            D(j) = randi(dmax,1); % dimension of this variable
            noise_id = randi(length(noise),1);
            if(noise_id==1)
                X{j} = randn(T,D(j));
            end
            if(noise_id==2)
                X{j} = (rand(T,D(j))-0.5);
            end
        else
            PA = find(G(:,j)==1);
            nPA = length(PA); % number of parents
            D(j) = randi(dmax,1);
            func_id = randi(length(func),1);
            noise_id = randi(length(noise),1);
            Xpa = [];
            for ii = 1:nPA
                Xpa = [Xpa X{PA(ii)}];
            end
            nPA2 = size(Xpa,2);
            if(func_id==1)
                X{j} = Xpa*ones(nPA2,D(j))*1.7/(nPA2+1);
            end
            if(func_id==2)
                X{j} = sin(Xpa*ones(nPA2,D(j)));
            end
            if(func_id==3)
                X{j} = cos(Xpa*ones(nPA2,D(j)));
            end
            if(func_id==4)
                X{j} = tanh(Xpa*ones(nPA2,D(j)));
            end
            if(noise_id==1)
                X{j} = X{j} + 0.4*randn(T,D(j));
            end
            if(noise_id==2)
                X{j} = X{j} + 0.5*(rand(T,D(j))-0.5);
            end
        end
    end
    
    XX = [];
    d_label = [];
    tmp = 0;
    for ii = 1:N
        XX = [XX X{ii}];
        d_label{ii} = [tmp+1:tmp+size(X{ii},2)];
        tmp = tmp+size(X{ii},2);
    end
    
    
    C = corr(XX); % roughly check the faithfulness by estimating the correlation matrix
    sign = 0;
    for ii = 1:N
        for jj = 1:N
            if(G(ii,jj)==1)
                tmpC = C(d_label{ii},d_label{jj});
                if(sum(abs(tmpC(:)))<0.4 | ~isempty(find(abs(tmpC(:))>0.9)))
                    sign = 1;
                    break;
                end
            end
        end
        if(sign)
            break;
        end
    end
    if(~sign)
        count = count+1;
        G_save{count} = G;
        Data_save{count} = XX;
        d_label_save{count} = d_label;
    end
end

save data_generate_multi G_save Data_save d_label_save

