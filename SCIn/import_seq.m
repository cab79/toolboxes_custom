function import_seq

orig_file = uigetfile('Sequence*.mat','Select file to take sequence from');
load(orig_file)
orig_seq=seq;

new_file = uigetfile('Sequence*.mat','Select file to add sequence to');
load(new_file)
try
    adapt=seq.adapttype;
end
seq=orig_seq;
try
    seq.adapttype=adapt;
end
save(new_file,'filename','seq','settings')