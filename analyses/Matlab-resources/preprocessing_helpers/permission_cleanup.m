function permission_cleanup(filepath)

command = ['chgrp -hR nbp ' fullfile(filepath)];
system(command);
command = ['chmod 771 -R ' fullfile(filepath)];
system(command);

end



