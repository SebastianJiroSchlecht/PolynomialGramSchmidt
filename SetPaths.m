% Add the pinned toolbox and supplied benchmark implementations.
projectDir = fileparts(mfilename('fullpath'));
assert(isfolder(fullfile(projectDir,'PEVDToolbox','General')),...
    'Initialize PEVDToolbox with git submodule update --init --recursive.');
addpath(projectDir);
addpath(fullfile(projectDir,'PEVDToolbox','General'),'-end');
addpath(fullfile(projectDir,'PEVDToolbox','Decompositions'),'-end');
addpath(fullfile(projectDir,'PQRD_code_john'),'-end');
addpath(fullfile(projectDir,'PQRD_code_faizan'),'-end');
addpath(fullfile(projectDir,'PQRD_code_fraser'),'-end');
