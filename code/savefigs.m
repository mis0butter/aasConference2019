function savefigs(alpha)
% Inputs 
%   alpha :     string 'alpha0' or 'alphaNot0' 

FolderName = fullfile('figures', alpha);   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Name'));
  set(0, 'CurrentFigure', FigHandle);
  saveas(FigHandle, fullfile(FolderName, FigName), 'pdf')
  saveas(FigHandle, fullfile(FolderName, FigName), 'eps')
  saveas(FigHandle, fullfile(FolderName, FigName), 'png')
%   saveas(FigHandle, fullfile(FolderName, FigName), 'fig')
end

end 