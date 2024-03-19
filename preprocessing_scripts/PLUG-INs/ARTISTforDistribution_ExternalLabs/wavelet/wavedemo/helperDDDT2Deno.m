function y = helperDDDT2Deno(meth,x,thr)
% This function helperDDDT2Deno is only in support of DualtreeExample.
% It may change in a future release.
%
%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 09-Nov-2012.
%   Last Revision: 24-Jul-2013.
%   Copyright 1995-2013 The MathWorks, Inc.
 

Depth = 2;
switch meth
    case 'dwt'
        Df = dtfilters('farras');
        dt = dddtree2('dwt',x,Depth,Df);
        cfs = dt.cfs;
        % loop through scales:
        for j = 1:Depth
            % loop through subbands
            for s = 1:3
                cfs{j}(:,:,s) = wthresh(cfs{j}(:,:,s),'s',thr);
            end
        end
        
    case 'realdt'
        FDf = dtfilters('FSfarras');
        Df  = dtfilters('qshift10');
        dt = dddtree2('realdt',x,Depth,FDf,Df);
        cfs = dt.cfs;
        % loop through scales:
        for j = 1:Depth
            % loop through subbands
            for s1 = 1:2
                for s2 = 1:3
                    cfs{j}(:,:,s2,s1) = wthresh(cfs{j}(:,:,s2,s1),'s',thr);
                end
            end
        end
        
    case 'cplxdt'
        FDf = dtfilters('FSfarras');
        Df   = dtfilters('qshift10');
        dt = dddtree2('cplxdt',x,Depth,FDf,Df);
        cfs = dt.cfs;
        % loop through scales:
        for j = 1:Depth
            % loop through subbands
            for s1 = 1:2
                for s2 = 1:3
                    C = cfs{j}(:,:,s2,s1,1) + 1i*cfs{j}(:,:,s2,s1,2);
                    C = wthresh(C,'s',thr);
                    cfs{j}(:,:,s2,s1,1) = real(C);
                    cfs{j}(:,:,s2,s1,2) = imag(C);
                end
            end
        end

    case 'ddt'
        Df = dtfilters('filters1');
        dt = dddtree2('ddt',x,Depth,Df);
        cfs = dt.cfs;
        % loop through scales
        for j = 1:Depth
            % loop through subbands
            for s = 1:8
                cfs{j}(:,:,s) = wthresh(cfs{j}(:,:,s),'s',thr);
            end
        end
        
    case 'realdddt'
        FDf = dtfilters('FSdoubledualfilt');
        Df = dtfilters('doubledualfilt');
        dt = dddtree2('realdddt',x,Depth,FDf,Df);
        cfs = dt.cfs;
        % loop through scales:
        for j = 1:Depth
            % loop through subbands
            for s1 = 1:2
                for s2 = 1:8
                    cfs{j}(:,:,s2,s1) = wthresh(cfs{j}(:,:,s2,s1),'s',thr);
                end
            end
        end
        
    case 'cplxdddt'
        FDf = dtfilters('FSdoubledualfilt');
        Df = dtfilters('doubledualfilt');        
        dt = dddtree2('cplxdddt',x,Depth,FDf,Df);
        cfs = dt.cfs;
        % loop through scales:
        for j = 1:Depth
            % loop through subbands
            for s1 = 1:2
                for s2 = 1:8
                    C = cfs{j}(:,:,s2,s1,1) + 1i*cfs{j}(:,:,s2,s1,2);
                    C = wthresh(C,'s',thr);
                    cfs{j}(:,:,s2,s1,1) = real(C);
                    cfs{j}(:,:,s2,s1,2) = imag(C);
                end
            end
        end
end
dt.cfs = cfs;
y = idddtree2(dt);

