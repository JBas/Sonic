sim: sim.c
	gcc sim.c `pkg-config --libs --cflags raylib` -o sim

save: sim_save.c
	gcc sim_save.c `pkg-config --libs --cflags raylib` -o save

clean:
	rm sim save
