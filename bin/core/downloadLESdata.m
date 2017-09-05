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
        error('Could not download LES data from URL. Check your internet settings.');
    end
end
% disp('Found the necessary LES files.');
end
