function check_folder(folder_name)
%CHECK_FOLDER Summary of this function goes here
%   Detailed explanation goes here
    current_folder = pwd();
    cd(exportpath())
    if not(isfolder(folder_name))
        mkdir(folder_name)
    end
    cd(current_folder)
end
