
%%% cut2wav.m

function [TDTdata] = cut2wav()

uiwait(msgbox('Select TDT data stored from FF2matsingle()', 'Instructions', "modal")); 

uiopen();

fs = TDTdata.streams.x465N.fs;

% Check if file has been already cut or was already registered only for the
% 5 min in the FC-box:

if round(length(TDTdata.streams.x465N.data)./fs) == 300 || round(length(TDTdata.streams.x465N.data)./fs) < 301
    mes = 'Data is already 300 sec long or is not the right file';
    waitfor(warndlg(mes, 'Please check your data'))

else

% Recortar señal de FF y modificar vectores en results:


idx_time0 = find(TDTdata.streams.Wav1.data > 3, 1, 'first');

fiveminGCAMP = TDTdata.streams.x465N.data(idx_time0:(idx_time0 + round(300*fs)));
fiveminIsos = TDTdata.streams.x405N.data(idx_time0:(idx_time0 + round(300*fs)));


TDTdata.streams.x465N.data = fiveminGCAMP;
TDTdata.streams.x405N.data = fiveminIsos;

uisave('TDTdata', 'fpdata.mat');


end % del condicional

end % del bucle