#

LN      =       ln -sf
MAKE    =       make -i -r
RM      =       rm -f

MODULES =  \
	module_stringutil.o \
	misc_definitions_module.o \
	parallel_module.o \
	hash_module.o \
	smooth_module.o \
	module_debug.o \
	bitarray_module.o \
	queue_module.o \
	list_module.o \
	gridinfo_module.o \
	interp_module.o \
	llxy_module.o \
	source_data_module.o \
	proc_point_module.o \
	process_tile_module.o \
	module_geogrid.o \
	module_real.o \
	datatype_module.o \
	module_metgrid.o \
	interp_option_module.o \
	process_domain_module.o \
	met_data_module.o \
	minheap_module.o \
	module_mergesort.o \
	module_mergesort.o \
	rotate_winds_module.o \
	storage_module.o \
	module_camdomain.o \
	module_geodomain.o \
	module_bdy_prep.o

OBJS    = cio.o read_geogrid.o

LIBPATHS =

LIBTARGET    =  prepare
TARGETDIR    =  ./
$(LIBTARGET) :  $(MODULES) $(OBJS)
		$(AR) $(ARFLAGS) ../main/$(LIBWRFLIB) $(MODULES) $(OBJS)

include ../configure.wrf

.c.o:
	$(RM) $@
	$(CC) -o $@ -c -D_UNDERSCORE $(CFLAGS) $*.c

clean:
	@ echo 'use the clean script'

# DEPENDENCIES : only dependencies after this line (don't remove the word DEPENDENCIES)


# this already built above :../dyn_em/module_initialize.real.o 

parallel_module.o: \
	../frame/module_dm.o \
	../frame/module_atm_communicator.o

cio.o:

read_geogrid.o:

module_debug.o: \
	cio.o \
	parallel_module.o

bitarray_module.o: \
	module_debug.o

list_module.o: \
	module_debug.o

queue_module.o: \
	module_debug.o

smooth_module.o: \
	parallel_module.o

gridinfo_module.o: \
	../share/module_llxy.o \
	misc_definitions_module.o \
	module_debug.o \
	../frame/module_driver_constants.o

interp_module.o: \
	bitarray_module.o \
	misc_definitions_module.o \
	module_debug.o \
	queue_module.o \
        ../frame/module_dm.o

llxy_module.o: \
	gridinfo_module.o \
	list_module.o \
	../share/module_llxy.o \
	module_debug.o \
	misc_definitions_module.o

source_data_module.o: \
	hash_module.o \
	read_geogrid.o \
	list_module.o \
	module_debug.o \
	misc_definitions_module.o \
	gridinfo_module.o \
	module_stringutil.o \
	llxy_module.o \
	../share/module_llxy.o

proc_point_module.o: \
	bitarray_module.o \
	hash_module.o \
	misc_definitions_module.o \
	module_debug.o \
	source_data_module.o \
	llxy_module.o \
	interp_module.o

process_tile_module.o: \
	bitarray_module.o \
	interp_module.o \
	llxy_module.o \
	smooth_module.o \
	misc_definitions_module.o \
	proc_point_module.o \
	queue_module.o \
	source_data_module.o \
	hash_module.o \
	../share/module_llxy.o \
	gridinfo_module.o

module_geogrid.o: \
	../share/module_llxy.o \
	gridinfo_module.o \
	llxy_module.o \
	list_module.o \
	module_debug.o \
	process_tile_module.o \
	source_data_module.o \
	parallel_module.o

module_real.o: \
	../frame/module_domain.o \
	../frame/module_configure.o \
	../dyn_em/module_big_step_utilities_em.o \
	../dyn_em/module_initialize_real.o \
	../share/module_optional_input.o \
	../share/module_bc.o \
	module_bdy_prep.o

module_bdy_prep.o: \
		../frame/module_domain.o \
		../frame/module_configure.o  \
		../frame/module_state_description.o \
		../frame/module_timing.o \
		../frame/module_dm.o \
		../dyn_em/module_initialize_real.o \
		../share/module_model_constants.o \
		../share/module_bc.o \
		../share/module_optional_input.o \
		../prep/parallel_module.o

datatype_module.o: bitarray_module.o module_debug.o

module_metgrid.o: \
	interp_option_module.o \
	module_debug.o \
	process_domain_module.o \
	../frame/module_domain.o \
	../frame/module_configure.o \
	module_camdomain.o

interp_option_module.o: list_module.o \
	misc_definitions_module.o \
	module_debug.o \
	module_stringutil.o \
	../share/module_llxy.o

met_data_module.o:

minheap_module.o: datatype_module.o

module_mergesort.o:

process_domain_module.o: bitarray_module.o \
	interp_module.o \
	interp_option_module.o \
	met_data_module.o \
	list_module.o \
	llxy_module.o \
	misc_definitions_module.o \
	module_debug.o \
	module_mergesort.o \
	parallel_module.o \
	rotate_winds_module.o \
	storage_module.o \
	../frame/module_domain.o \
	../frame/module_configure.o \
	module_camdomain.o \
        queue_module.o

rotate_winds_module.o: bitarray_module.o \
	llxy_module.o \
	misc_definitions_module.o \
	module_debug.o \
	../frame/module_driver_constants.o

storage_module.o: datatype_module.o \
	minheap_module.o \
	misc_definitions_module.o \
	module_debug.o \
	parallel_module.o \
	module_stringutil.o

module_camdomain.o: parallel_module.o \
	../frame/module_dm.o \
	../frame/module_domain.o

module_geodomain.o: 

# DO NOT DELETE
