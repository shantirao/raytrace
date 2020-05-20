function result= testSurfaceLocal
%% test surfaceLocal

result = true;
result = result | doTest([0 0 1], [1 0 0; 0 1 0]);
result = result | doTest([1 0 0], [0 0 -1; 0 1 0]);
result = result | doTest([0 1 0], [1 0 0; 0 0 -1]);

%% special cases

%%
function s = doTest(d,r)
    t = surfaceLocal(struct('direction',d));
    s = isequal(r,t);
    if ~s
        disp(['error in: ' sprintf('%d ',d) ' out: ' sprintf('%d ',t)  ' expected: ' sprintf('%d ',r)]);
    end
end
end