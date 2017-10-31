function [ fn1,fn2 ] = downloadLESdata( WFSimFolder, meshingName )
fn1 = [WFSimFolder 'data_LES/' meshingName '_meshing.mat'];
fn2 = [WFSimFolder 'data_LES/' meshingName '_data.mat'];
if exist(fn1) ~= 2 | exist(fn2) ~= 2
    disp('File not found. Attempting to downloading from URL..');
    try
        url1 = ['http://homepage.tudelft.nl/h5h68/data_LES/' meshingName '_meshing.mat'];
        url2 = ['http://homepage.tudelft.nl/h5h68/data_LES/' meshingName '_data.mat'];
        outFn1 = websave(fn1,url1);
        outFn2 = websave(fn2,url2);
    catch
        delete([fn1 '*']); delete([fn2 '*']); % Delete junk files from download
        disp(' ');
        error(sprintf(['Could not download LES data from URL. Check your internet settings.\n' ...
            'Also, make sure your MATLAB version supports the function "websave(..)".\n\n' ...
            'If this is not the case, please download the following files manually\n' ...
            'And place them in the data_LES folder:\n\n' ...
            'File 1: ' url1 '\n' ...
            'File 2: ' url2 '\n\n' ...
            'Reach out to us on Github or through email if the error persists.\n']))
    end
end
% disp('Found the necessary LES files.');
end
