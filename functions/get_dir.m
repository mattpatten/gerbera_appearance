function [varargout] = get_dir(imgAbbr, varargin)

% This function is a way to put the entire directory structure in one file, such that it need not be mentioned
% in code anywhere else. It takes a variable number of inputs and outputs so you specify which directories you 
% want (preceded by the dataset abbreviation), and its output variable name.
%
% Inputs:
%    imgAbbr     - Abbreviated description of the dataset of images to use.
%    <varargin>  - A string, or set of strings containing any directories you want to retrieve.
%                  Strings should be the text preceding 'Dir' in variable names below.
%                  Currently: root, output, origStim, procStim, img, mask, proc, rem, data
%
% Outputs:
%    <varargout> - The variable names of any directories you are retrieving.
%
% e.g.1,
%   root = get_dir([],'root'); %gets the root directory
%
% e.g.2, 
%   [imgDir, maskDir] = get_dir('google','img','mask')
%   %retrieves image and mask directories from the 'google' dataset
%
% Created by Matt Patten in Nov, 2018


%default value if imgAbbr is not specified
if ~exist('imgAbbr','var'), imgAbbr = ''; end

%define folders (that don't require imgAbbr)
rootDir     = [fileparts(mfilename('fullpath')) filesep '..' filesep]; %base directory for experiment
dataDir     = [rootDir 'data'    filesep];                             %base directory of data files
origStimDir = [rootDir 'stimuli' filesep 'Original'  filesep];         %directory where original (untouched) image dataset is stored 
procStimDir = [rootDir 'stimuli' filesep 'Processed' filesep];         %directory where different sets of output stimulus images are stored
imgPropDir  = [rootDir 'img_properties'  filesep];                     %directory containing things related to gathering the image properties
databaseDir = '.\..\PostgreSQL\pg_output\'; 

%define folders that do require the imgAbbr parameter
outputDir  = [rootDir     'output'  filesep imgAbbr  filesep];         %base for outputs
imgDir     = [procStimDir 'images_' imgAbbr filesep];                  %directory where individual images are located
maskDir    = [imgDir      'masks'           filesep];                  %directory containing images of masks 
procDir    = [imgDir      'proc'            filesep];                  %directory containing images that display generation of mask
remDir     = [imgDir      'rem'             filesep];                  %directory containing images that were removed from the analysis
remMaskDir = [remDir      'masks'           filesep];                  %directory containing masks that were removed from the analysis

%return the directory/directories that were asked for inputs should be strings, with whatever is before Dir in the above 
%variable names (e.g., 'stim' for stimDir)
for i=1:nargin-1
    varargout{i} = eval([varargin{i} 'Dir']);
end

end