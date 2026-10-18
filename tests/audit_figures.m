function audit_figures
% Run each paper plotting script in a fresh workspace, retaining all failures.
SetPaths;
set(groot,'defaultFigureVisible','off');
scripts = {'TSP26_PolyGramSchmidt_Figs12','TSP26_PolyGramSchmidt_Figs34',...
    'TSP26_PolyGramSchmidt_Figs56','TSP26_PolyGramSchmidt_Fig7'};
diary('reproduction-log.txt');
cleanup = onCleanup(@() diary('off'));
failures = {};
for k = 1:numel(scripts)
    fprintf('\nRUN %s\n',scripts{k});
    try
        run_one(scripts{k});
        figures = findall(groot,'Type','figure');
        for n = 1:numel(figures)
            print(figures(n),sprintf('audit-%s-%d.png',scripts{k},figures(n).Number),'-dpng','-r100');
        end
        close all;
        fprintf('PASS %s\n',scripts{k});
    catch exception
        fprintf('FAIL %s\n%s\n',scripts{k},getReport(exception));
        failures{end+1} = scripts{k}; %#ok<AGROW>
    end
end
assert(isempty(failures),'Figure reproduction failed: %s',strjoin(failures,', '));
end

function run_one(scriptName)
run(scriptName);
end
