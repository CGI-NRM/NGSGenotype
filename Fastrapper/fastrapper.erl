%%% Fastrapper: The FASTQ Bootstrapper
-module(fastrapper).
-export([reader_subworker/2, reader_worker/3, supervisor/3, random_sampler/3, start/1]).

%%% Functions for reading:
file_reader(File_name) ->
	{ok, File} = file:read_file(File_name),
	Inflated_data = zlib:gunzip(File),
	binary:split(Inflated_data, [<<"\n">>], [global]).

tuplify([], Output) ->
	list_to_tuple(Output);

tuplify([_], Output) ->
	list_to_tuple(Output);

tuplify([Row1, Row2, Row3, Row4|T], Output) ->
	tuplify(T, [{Row1, Row2, Row3, Row4}|Output]).

reader_subworker(Path, PID) ->
	File = file_reader(Path),
	PID ! {fileDone, tuplify(File, [])}.

reader_worker({Path, Name}, _, PID) ->
	spawn(fastrapper, reader_subworker, [Path, self()]), % Read forward file in the background
	[Short_path|_] = string:replace(Path, "_1.fq.gz", ""),
	[Short_name|_] = string:replace(Name, "_1.fq.gz", ""),
	File_R = file_reader(Short_path ++ "_2.fq.gz"),
	Tuplist_R = tuplify(File_R, []),
	receive
		{fileDone, Tuplist_F} ->
			ok
	end,
	PID ! {part, {Short_name, Tuplist_F, Tuplist_R}}.

mass_reader(File_folder) ->
	{ok, File_names} = file:list_dir(File_folder),
	Fq_tester = fun(X) -> (re:run(X, "_1.fq.gz", [{capture, none}]) == match) end,
	Fq_names = lists:filter(Fq_tester, File_names),
	Fq_paths = lists:map(fun(X) -> File_folder ++ "/" ++ X end, Fq_names),
	io:fwrite("Forward files found:~n"),
	lists:foreach(fun(X) -> io:fwrite("~s~n", [X]) end, Fq_paths),
	
	PID = spawn(fastrapper, supervisor, [[], length(Fq_paths), self()]),
	worker_launcher(reader_worker, lists:zip(Fq_paths, Fq_names), 0, PID),
	receive
		{allParts, Package} ->
			Package
	end.

%%% Functions for sampling:
random_sampler({Name, Pop_for, Pop_rev}, Size, PID) ->
	File_name = Name ++ "_" ++ integer_to_list(Size) ++ "_bootstraps",
	Length = tuple_size(Pop_for),
	Draws = [rand:uniform(Length) || _ <- lists:seq(1, Size)],
	PID ! {part, {File_name, [{element(X, Pop_for), element(X, Pop_rev)} || X <- Draws]}}.

mass_sampler(File_list, Sample_size) ->
	PID = spawn(fastrapper, supervisor, [[], length(File_list), self()]),
	worker_launcher(random_sampler, File_list, Sample_size, PID),

	receive
		{allParts, Package} ->
			Package
	end.

%%% Functions for writing:
f_r_file([], File_F, File_R) ->
	{File_F, File_R};

f_r_file([{{F1, F2, F3, F4}, {R1, R2, R3, R4}}|T], File_F, File_R) ->
	File_F_Out = io_lib:fwrite("~s~n~s~n~s~n~s~n", [F1, F2, F3, F4]) ++ File_F,
	File_R_Out = io_lib:fwrite("~s~n~s~n~s~n~s~n", [R1, R2, R3, R4]) ++ File_R,
	f_r_file(T, File_F_Out, File_R_Out).

gzipper(Name, File) ->
	Gz_file = zlib:gzip(File),
	file:write_file(Name ++ ".gz", Gz_file).

mass_writer([]) ->
	allWritten;

mass_writer([{Name, Samplings}|T]) ->
	Forward = Name ++ "_1.fq",
	Reverse = Name ++ "_2.fq",
	{File_F, File_R} = f_r_file(Samplings, [], []),
	gzipper(Forward, File_F),
	gzipper(Reverse, File_R),
	mass_writer(T).

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

%%% Start program:
start(Args) ->
	[File_folder, Bs_string] = Args,
	Bootstraps = list_to_integer(Bs_string),
	File_list = mass_reader(File_folder),
	io:fwrite("~nFiles loaded.~n"),
	Sample_list = mass_sampler(File_list, Bootstraps),
	io:fwrite("Files bootstrapped.~n"),
	mass_writer(Sample_list),
	io:fwrite("Files written.~n").
