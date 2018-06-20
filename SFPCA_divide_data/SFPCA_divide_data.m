  subjects = [1,2,3,4,5,6,8,9,10,11,12,15,16,17,18,19,21];
  tic;
  load('scanidx')
  for s=1:length(subjects)
    subject = subjects(s);
    filename = ['fmridata',num2str(subject),'cnew'];
    load(filename)
%     filename = ['scanidx',num2str(subject)];
%     load(filename)
    for q=1:81
      temp_q = temp(:,:,:,scan_idx(q,:)); 
      filename = ['voxels_',num2str(subject),'_',num2str(q)];
      save(filename, 'temp_q')
    end
  end
  toc;