function [U, V, lambdas, iterLog] = mMultiNMF_l21p1_ard(X, W, D, L, Wt, Dt, Lt, initU, initV, options)
% with ARD regularization on U and V

%viewNum = length(X);

Rounds = options.rounds;
maxInIter = options.maxIter;
minInIter = options.minIter;

U = initU;
V = initV;

clear initU;
clear initV;

%
if isfield(options,'ardUsed') && options.ardUsed>0
    disp('mMultiNMF_l21p1 with ard...');
    %hyperLam = zeros(viewNum,1);
    %lambdas = cell(viewNum,1);
    eta = options.eta;
    %for vi=1:viewNum
    [mFea,nSmp] = size(X);
    lambdas = sum(U) / mFea;
    
    hyperLam = eta * sum(sum(X.^2)) / (mFea*nSmp*2);
    %end
else
    disp('mMultiNMF_l21p1...');
end

oldL = Inf;
j = 0;
iterLog = 0;
restartJ = 0; % added by hmli
while j < Rounds
    % calculate current objective function value
    j = j + 1;
    
    tmpl21 = zeros(size(V));
    L1 = 0;
    ardU = 0;
    tmp1 = 0;
    tmp2 = 0;
    tmp3 = 0;
    
    % for i = 1:viewNum
    [mFea,nSmp] = size(X);
    
    tmpl21 = tmpl21 + V.^2;
    
    if isfield(options,'alphaS1')
        tmpNorm1 = sum(V,1);
        tmpNorm2 = sqrt(sum(V.^2,1));
        L1 = L1 + options.alphaS1 * sum(tmpNorm1./max(tmpNorm2,eps));
    end
    
    % ard term for U
    if isfield(options,'ardUsed') && options.ardUsed>0
        su = sum(U);
        su(su==0) = 1;
        ardU = ardU + sum(log(su))*mFea*hyperLam;
    end
    
    tmpDf = (X-U*V').^2;
    tmp1 = tmp1 + sum(tmpDf(:));
    
    if isfield(options,'alphaL')
        dVi = double(V');
        tmp2 = tmp2 + dVi * L .* dVi;
    end
    
    if isfield(options,'alphaLT')
        dUi = double(U');
        tmp3 = tmp3 + dUi * Lt .* dUi;
    end
    
    L21 = options.alphaS21 * sum(sum(sqrt(tmpl21))./max(sqrt(sum(tmpl21)),eps));
    Ldf = tmp1;
    Lsl = sum(tmp2(:));
    Ltl = sum(tmp3(:));
    
    logL = L21 + ardU + Ldf + Lsl + Ltl + L1;
    
    iterLog(end+1) = logL;
    %disp(['  round:',num2str(j),' logL:',num2str(logL)]);
    %     disp(['  round:',num2str(j),' logL:',num2str(logL),',dataFit:',num2str(Ldf)...
    %           ',spaLap:',num2str(Lsl),',temLap:',num2str(Ltl),',L21:',num2str(L21),...
    %           ',L1:',num2str(L1),',ardU:',num2str(ardU)]);
    
    %    if (oldL < logL)
    %        % modified by hmli
    %         if restartJ == 0
    %             restartJ = j;
    %             consRestart = 1;
    %         elseif restartJ == j
    %             consRestart = consRestart + 1;
    %             if consRestart>=3
    % 				U = oldU;
    % 				V = oldV;
    % 				logL = oldL;
    %
    %                 break;
    %             end
    %         else
    %             restartJ = j;
    %             consRestart = 1;
    %         end
    %
    %         j = j - 1;
    %         %disp('restrart this iteration');
    %    else
    if j>5 && (oldL-logL)/max(oldL,eps)<options.error
        break;
    end
    %    end
    
    oldU = U;
    oldV = V;
    oldL = logL;
    
    %for i=1:viewNum%%paralell
    [mFea,nSmp] = size(X);
    
    iter = 0;
    oldInLogL = inf;
    
    fixl2 = zeros(size(V));
    %for vi = 1:viewNum
    %if vi~=i
    %fixl2 = fixl2 + V.^2;
    %end
    %end
    
    while iter<maxInIter
        iter = iter + 1;
        
        % ===================== update V ========================
        XU = X'*U;
        UU = U'*U;
        VUU = V*UU;
        
        tmpl2 = fixl2 + V.^2;
        if options.alphaS21>0
            tmpl21 = sqrt(tmpl2);
            tmpl22 = repmat(sqrt(sum(tmpl2,1)),nSmp,1);
            tmpl21s = repmat(sum(tmpl21,1),nSmp,1);
            posTerm = V ./ max(tmpl21.*tmpl22,eps);
            negTerm = V .* tmpl21s ./ max(tmpl22.^3,eps);
            
            VUU = VUU + 0.5 * options.alphaS21 * posTerm;
            XU = XU + 0.5 * options.alphaS21 * negTerm;
        end
        
        if isfield(options,'alphaL')
            WV = W * double(V);
            DV = D * double(V);
            
            XU = XU + WV;
            VUU = VUU + DV;
        end
        
        if isfield(options,'alphaS1')
            sV = max(repmat(sum(V),nSmp,1),eps);
            normV = sqrt(sum(V.^2));
            normVmat = repmat(normV,nSmp,1);
            posTerm = 1./max(normVmat,eps);
            negTerm = V.*sV./max(normVmat.^3,eps);
            
            XU = XU + 0.5*options.alphaS1*negTerm;
            VUU = VUU + 0.5*options.alphaS1*posTerm;
        end
        
        V = V.*(XU./max(VUU,eps));
        
        prunInd = sum(V~=0)==1;
        if any(prunInd)
            V(:,prunInd) = zeros(nSmp,sum(prunInd));
            U(:,prunInd) = zeros(mFea,sum(prunInd));
        end
        
        % ==== normalize U and V ====
        [U,V] = Normalize(U, V);
        
        % ===================== update U =========================
        XV = X*V;
        VV = V'*V;
        UVV = U*VV;
        
        if isfield(options,'ardUsed') && options.ardUsed>0 % ard term for U
            posTerm = 1./max(repmat(lambdas,mFea,1),eps);
            UVV = UVV + posTerm*hyperLam;
        end
        
        if isfield(options,'alphaLT')
            WU = Wt * double(U);
            DU = Dt * double(U);
            
            XV = XV + WU;
            UVV = UVV + DU;
        end
        
        U = U.*(XV./max(UVV,eps));
        %U{i}(U{i}<1e-6) = 0;
        
        prunInd = sum(U)==0;
        if any(prunInd)
            V(:,prunInd) = zeros(nSmp,sum(prunInd));
            U(:,prunInd) = zeros(mFea,sum(prunInd));
        end
        
        % update lambda
        if isfield(options,'ardUsed') && options.ardUsed>0
            lambdas = sum(U) / mFea;
        end
        % ==== calculate partial objective function value ====
        inTl = 0;
        inSl = 0;
        LardU = 0;
        LL1 = 0;
        
        inDf = (X-U*V').^2;
        
        if isfield(options,'alphaLT')
            dUi = double(U');
            inTl = dUi * Lt .* dUi;
        end
        if isfield(options,'alphaL')
            dVi = double(V');
            inSl = dVi * L .* dVi;
        end
        if isfield(options,'ardUsed') && options.ardUsed>0
            % ard term for U
            su = sum(U);
            su(su==0) = 1;
            LardU = sum(log(su))*mFea*hyperLam;
        end
        inL21 = zeros(size(V));
        if options.alphaS21>0
            %for vi=1:viewNum
                inL21 = inL21 + V.^2;
            %end
        end
        if isfield(options,'alphaS1')
            tmpNorm1 = sum(V,1);
            tmpNorm2 = sqrt(sum(V.^2,1));
            LL1 = options.alphaS1 * sum(tmpNorm1./max(tmpNorm2,eps));
        end
        
        inL21 = sum(sqrt(inL21))./max(sqrt(sum(inL21)),eps);
        LDf = sum(inDf(:));
        LTl = sum(inTl(:));
        LSl = sum(inSl(:));
        LL21 = options.alphaS21 * sum(inL21(:));
        
        inLogL = LDf + LTl + LSl + LardU + LL21 + LL1;
        
        if iter>minInIter && abs(oldInLogL-inLogL)/max(oldInLogL,eps)<options.error
            break;
        end
        oldInLogL = inLogL;
    end
end
end % function


function [U, V] = Normalize(U, V)
[U,V] = NormalizeUV(U, V, 1, 1);
end


function [U, V] = NormalizeUV(U, V, NormV, Norm)
nSmp = size(V,1);
mFea = size(U,1);
if Norm == 2
    if NormV
        norms = sqrt(sum(V.^2,1));
        norms = max(norms,eps);
        V = V./repmat(norms,nSmp,1);
        U = U.*repmat(norms,mFea,1);
    else
        norms = sqrt(sum(U.^2,1));
        norms = max(norms,eps);
        U = U./repmat(norms,mFea,1);
        V = V.*repmat(norms,nSmp,1);
    end
else
    if NormV
        %norms = sum(abs(V),1);
        norms = max(V);
        norms = max(norms,eps);
        V = V./repmat(norms,nSmp,1);
        U = U.*repmat(norms,mFea,1);
    else
        %norms = sum(abs(U),1);
        norms = max(U);
        norms = max(norms,eps);
        U = U./repmat(norms,mFea,1);
        V = bsxfun(@times, V, norms);
    end
end
end


