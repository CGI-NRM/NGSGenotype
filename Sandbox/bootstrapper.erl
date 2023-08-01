-module(bootstrapper).
-export([file_reader/1, sampling_supervisor/3, random_sampler/3, make_sample/2, test_sampling/0, test_reading/0, test_writing/0, test_readwrite/0, start/0]).

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

%%% Functions for sampling:
sampling_supervisor(Package, 0, PID) ->
	PID ! {allSamples, Package};

sampling_supervisor(Package, Chunks, PID) ->
	receive
		{part, Samples} ->
			sampling_supervisor([Samples|Package], Chunks - 1, PID)
	end.

random_sampler(Population, Size, PID) ->
	Length = length(Population),
	Pop_tuple = list_to_tuple(Population), % Since indexing tuples are way faster
	PID ! {part, [element(rand:uniform(Length), Pop_tuple) || _ <- lists:seq(1, Size)]}.

worker_launcher([], _, _) ->
	allSpawned;

worker_launcher([H|T], Sample_size, PID) ->
	spawn(bootstrapper, random_sampler, [H, Sample_size, PID]),
	worker_launcher(T, Sample_size, PID).

make_sample(File_list, Sample_size) ->
	PID = spawn(bootstrapper, sampling_supervisor, [[], length(File_list), self()]),
	worker_launcher(File_list, Sample_size, PID),

	receive
		{allSamples, Package} ->
			Package
	end.

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
	ok.
