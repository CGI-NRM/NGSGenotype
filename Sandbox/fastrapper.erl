%%% Fastrapper: The FASTQ Bootstrapper
-module(fastrapper).
-export([file_reader/1, reader_worker/3, supervisor/3, random_sampler/3, make_sample/2, test_sampling/0, test_reading/0, test_writing/0, test_readwrite/0, start/0]).

%%% Functions for reading and writing:
file_reader(File_name) ->
	{ok, File} = file:read_file(File_name),
	Inflated_data = zlib:gunzip(File),
	binary:split(Inflated_data, [<<"\n">>], [global]).

tuplify([], Output) ->
	Output;

tuplify([_], Output) ->
	Output;

tuplify([Row1, Row2, Row3, Row4|T], Output) ->
	tuplify(T, [{Row1, Row2, Row3, Row4}|Output]).

file_writer([], _) ->
	allWritten;

file_writer([{R1, R2, R3, R4}|T], Name) ->
	file:write_file(Name, io_lib:fwrite("~s~n~s~n~s~n~s~n", [R1, R2, R3, R4]), [append]),
	file_writer(T, Name).

gzipper(Name) ->
	{ok, File} = file:read_file(Name),
	Gz_file = zlib:gzip(File),
	file:write_file(Name ++ ".gz", Gz_file),
	file:delete(Name).

reader_worker(Path, _, PID) ->
	File_list = file_reader(Path),
	Tuplist = tuplify(File_list, []),
	PID ! {part, Tuplist}.	

mass_reader(File_folder) ->
	{ok, File_names} = file:list_dir(File_folder),
	Fq_tester = fun(X) -> (re:run(X, ".fq.gz", [{capture, none}]) == match) end,
	Fq_names = lists:filter(Fq_tester, File_names),
	Fq_paths = lists:map(fun(X) -> File_folder ++ "/" ++ X end, Fq_names),
	io:fwrite("~p~n", [Fq_paths]),
	
	PID = spawn(fastrapper, supervisor, [[], length(Fq_paths), self()]),
	worker_launcher(reader_worker, Fq_paths, 0, PID),
	receive
		{allParts, Package} ->
			Package
	end.

mass_writer([], _) ->
	allWritten;

mass_writer([File|File_tail], [Name|Name_tail]) ->
	file_writer(File, Name),
	gzipper(Name),
	mass_writer(File_tail, Name_tail).

%%% Functions for sampling:
random_sampler(Population, Size, PID) ->
	Length = length(Population),
	Pop_tuple = list_to_tuple(Population), % Since indexing tuples are way faster
	PID ! {part, [element(rand:uniform(Length), Pop_tuple) || _ <- lists:seq(1, Size)]}.

make_sample(File_list, Sample_size) ->
	PID = spawn(fastrapper, supervisor, [[], length(File_list), self()]),
	worker_launcher(random_sampler, File_list, Sample_size, PID),

	receive
		{allParts, Package} ->
			Package
	end.

%%% Functions for spawning:
supervisor(Package, 0, PID) ->
	PID ! {allParts, Package};

supervisor(Package, Chunks, PID) ->
	receive
		{part, Part} ->
			supervisor([Part|Package], Chunks - 1, PID)
	end.

worker_launcher(_, [], _, _) ->
	allSpawned;

worker_launcher(Function, [H|T], Setting, PID) ->
	spawn(fastrapper, Function, [H, Setting, PID]),
	worker_launcher(Function, T, Setting, PID).

%%% Test functions:
test_sampling() ->
	Task_list = [lists:seq(1, 1000000) || _ <- lists:seq(1, 6)],
	Results = make_sample(Task_list, 100000),
	io:fwrite("~p~n", [Results]).

test_reading() ->
	List = file_reader("../Raw_data/S646_142_EKDL230001504-1A_HNHKMDSX5_L3_1.fq.gz"),
	Tuplist = tuplify(List, []),
	Results = make_sample([Tuplist, Tuplist, Tuplist, Tuplist], 10000),
	io:fwrite("~p~n", [Results]).

test_writing() ->
	Fake_file = [{"@Potato1", "AACTGTCACG", "+", "FFFFFFFFFF"}, {"@Potato2", "AACTGTCACG", "+", "FFFFFFFFFF"}, {"@Potato3", "AACTGTCACG", "+", "FFFFFFFFFF"}],
	file_writer(Fake_file, "out.txt"),
	gzipper("out.txt").

test_readwrite() ->
	List = file_reader("../Raw_data/S646_142_EKDL230001504-1A_HNHKMDSX5_L3_1.fq.gz"),
	Tuplist = tuplify(List, []),
	Results = make_sample([Tuplist], 10000),
	file_writer(lists:nth(1, Results), "sample1.fq"),
	gzipper("sample1.fq").

%%% Start program:
start() ->
	File_folder = "../Raw_data",
	File_list = mass_reader(File_folder),
	Sample_list = make_sample(File_list, 10000),
	mass_writer(Sample_list, [integer_to_list(X) ++ ".txt" || X <- lists:seq(1, length(Sample_list))]),
	ok.
	% todo: read, sample and write all files with forward and reverse synched
	% todo: make the random sample file names correspond to the original
