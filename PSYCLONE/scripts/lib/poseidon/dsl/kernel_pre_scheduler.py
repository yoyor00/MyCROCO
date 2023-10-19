##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
from copy import deepcopy
# internal
from ..dsl.helper import *
from ..base.render_graph import RenderGraph
# psyclone
from psyclone.psyir.nodes import Node
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.transformations import TransformationError

##########################################################
class KernelPreScheduler:

    def __init__(self):
        self.files = {}
        self.routines = {}
        self.calls = {}
        self.meta = {}

    def add_file(self, fname: str, free_form=False):
        reader = FortranReader()
        psyir = reader.psyir_from_file(fname, free_form=free_form)
        self.add_psyir(fname, psyir)

    def add_psyir(self, fname: str, psyir: Node):
        self.files[fname] = {
            'psyir': psyir,
            'kernels': extract_kernels_from_psyir(psyir)
        }

    def make_loop_splitting_assignements(self):
        for fname, file in self.files.items():
            for kernel in file['kernels'].kernels:
                try:
                    kernel.make_loop_splitting_assignements()
                except TransformationError as error:
                    print(f"Fail to split kernel : {error}")
            file['kernels'] = extract_kernels_from_psyir(file['psyir'])
        #####
        #for kernel in file['kernels'].kernels:
        #    print(kernel.render_summary(orig=False, fused=False))
        #raise Exception("sdlfksdlkfj")

    def add_code(self, fname: str, code: str, free_form=False):
        reader = FortranReader()
        psyir = reader.psyir_from_source(code, free_form=free_form)
        self.add_psyir(fname, psyir)

    def save(self, fname: str):
        with open(fname, 'w+') as fp:
            json.dump(self.meta, fp, indent=True)

    def print(self):
        print(json.dumps(self.meta, indent=True))

    def load(self, fname: str):
        with open(fname, 'r') as fp:
            self.meta = json.load(fp)

    def print_codes(self):
        for fname, file in self.files.items():
            print(f"============ {fname} ==================")
            print(FortranWriter()(file['psyir']))

    def extract_routines(self):
        # extract routines
        self.routines = {}
        for fname, file in self.files.items():
            for routine in file['psyir'].walk(Routine):
                self.routines[routine.name] = routine

    def extract_calls(self):
        # extract routines
        self.calls = {}
        for fname, file in self.files.items():
            for call in file['psyir'].walk(Call):
                self.calls[call.routine.name] = call

    def apply(self):
        # extract
        self.extract_routines()
        self.extract_calls()
        self.assign_loop_ids()

        # apply meta
        for name, routine in self.routines.items():
            if name == 'step2d':
                continue
            routine.poseidon_previous_state = self.meta['routines'][name]['prev']
            routine.poseidon_end_state = self.meta['routines'][name]['end']

        # fill calls
        for name, call in self.calls.items():
            call.poseidon_previous_state = self.meta['calls'][name]['prev']
            call.poseidon_end_state = self.meta['calls'][name]['end']

        # apply loops infos
        for name, routine in self.routines.items():
            loops = routine.walk(Loop, stop_type=Loop)
            for loop in loops:
                id = loop.poseidon_func_loop_id
                print(f"{name}[{id}]")
                print(self.meta['loops'][name])
                loop.poseidon_acc_async_stream = self.meta['loops'][name][str(id)]['stream']
                loop.poseidon_acc_async_wait = self.meta['loops'][name][str(id)]['wait']

    def assign_loop_ids(self):
        # assign IDS to the loops
        for name, routine in self.routines.items():
            loops = routine.walk(Loop, stop_type=Loop)
            for id, loop in enumerate(loops):
                loop.poseidon_func_loop_id = id

    def get_child_kernels(self, routine: Routine):
        ordered_kernels = []
        for node in routine.walk((Loop, Call), stop_type=Loop):
            if isinstance(node, Loop):
                ordered_kernels.append(node.poseidon_kernel)
            elif isinstance(node, Call):
                ordered_kernels = ordered_kernels + self.get_child_kernels(self.routines[node.routine.name])
        return ordered_kernels

    def calculate2(self):
        # extract
        self.extract_routines()
        self.extract_calls()
        self.assign_loop_ids()

        # search roots
        self.roots = []
        for name, routine in self.routines.items():
            if not name in self.calls:
                self.roots.append(routine)

        # iterate over the routines
        #print(self.roots)
        root = self.roots[0]
        ordered_kernels = self.get_child_kernels(root)

        # loop over all
        streams = []
        var_stream = {}
        var_last_stream = {}

        # assign vars to streams
        for kernel in ordered_kernels:
            for var in kernel.out_vars:
                if not var in var_last_stream:
                    var_stream[var] = len(streams)
                    var_last_stream[var] = len(streams)
                    #if len(streams) < 8:
                    streams.append(len(streams))

        for kernel in ordered_kernels:
            kernel.schedule(True, var_stream, var_last_stream, streams)

        print(var_last_stream)
        old = var_last_stream.copy()

        # make second pass to handle multiple loops
        #for kernel in ordered_kernels:
        #    kernel.schedule(True, var_last_stream, streams)

        print(var_last_stream)
        for var in var_last_stream:
            if old[var] != var_last_stream[var]:
                print(f"diff : {{{var} = {old[var]} => {var_last_stream[var]}}}")
        assert var_last_stream == old

        # draw
        graph = RenderGraph()
        st = {}
        #for i in range(33):
        #    k = graph.add_loop(f"Stream {i}")
        #    st[i] = k
        for i, kernel in enumerate(ordered_kernels):
            out_vars = ','.join(kernel.make_uniq_sorted(kernel.out_vars))
            in_vars = ','.join(kernel.make_uniq_sorted(kernel.in_vars))

            k = graph.add_kernel(f"Kernel_{i}\nst={kernel.acc_async_stream}\nin={in_vars}\nout={out_vars}")
            if kernel.acc_async_stream in st:
                graph.add_link(st[kernel.acc_async_stream], k)
            st[kernel.acc_async_stream] = k
            if kernel.acc_async_wait != None:
                for dep in kernel.acc_async_wait:
                    if dep in st:
                        graph.add_link(st[dep], k)
        #graph.render_as_image("/home/svalat/tmp.png", "png")

        # store
        for name, routine in self.routines.items():
            loops = routine.walk(Loop, stop_type=Loop)
            if loops:
                id = 1
                routine.poseidon_previous_state = loops[0].poseidon_end_state
                while id <= len(loops):
                    try:
                        routine.poseidon_end_state = loops[-id].poseidon_end_state
                        break
                    except:
                        id += 1

        # gen summary
        meta = {
            'routines': {},
            'calls': {},
            'resync': [],
            'loops': {},
        }

        # build resync
        for var in var_stream:
            if var_last_stream[var] != var_stream[var]:
                meta['resync'].append({
                    'stream': var_stream[var], 
                    'wait': var_last_stream[var]
                })

        # fill routines
        cnt_streams = 0
        for name, routine in self.routines.items():
            print(f"Func : {name}()")
            if name == 'step2d':
                continue
            entry = {
                'prev': routine.poseidon_previous_state,
                'end': routine.poseidon_end_state
            }
            meta['routines'][name] = entry
            cnt_streams += len(routine.poseidon_end_state['streams'])

        # fill calls
        for name, call in self.calls.items():
            entry = {
                'prev': self.routines[name].poseidon_previous_state,
                'end': self.routines[name].poseidon_end_state
            }
            meta['calls'][name] = entry

        # fill loop ids
        for name, routine in self.routines.items():
            meta['loops'][name] = {}
            loops = routine.walk(Loop, stop_type=Loop)
            for loop in loops:
                wait = loop.poseidon_kernel.acc_async_wait
                if wait != None:
                    wait = list(wait)
                meta['loops'][name][str(loop.poseidon_func_loop_id)] = {
                    "stream": loop.poseidon_kernel.acc_async_stream,
                    "wait": wait,
                    #"check": FortranWriter()(loop)
                }

        # final result
        self.meta = meta

    def calculate(self, recurse=10):
        # extract
        self.extract_routines()
        self.extract_calls()

        # start converging
        old_meta = None
        meta = {}
        cnt = 0
        while meta != old_meta:
            # copy to compare on next iteration
            old_meta = deepcopy(meta)

            # calc async deps
            for fname, file in self.files.items():
                file['kernels'].calc_acc_async_wait(True)

            # apply previous loop on call & copy to routine entry
            for fname, file in self.files.items():
                for call in file['psyir'].walk(Call):
                    loop = get_previous(call, Loop)
                    if loop != None:
                        call.poseidon_previous_state = loop.poseidon_end_state
                        self.routines[call.routine.name].poseidon_previous_state = call.poseidon_previous_state

            # apply routine end on call
            for fname, file in self.files.items():
                for call in file['psyir'].walk(Call):
                    call.poseidon_end_state = self.routines[call.routine.name].poseidon_end_state.copy()
                    self.calls[call.routine.name] = call

            # apply unknown
            for fname, file in self.files.items():
                for call in file['psyir'].walk(Call):
                    loop = get_previous(call, Loop)
                    if loop == None:
                        call.poseidon_previous_state = call.poseidon_end_state

            # self loop missings
            for name, routine in self.routines.items():
                try:
                    _ = routine.poseidon_previous_state
                except:
                    try:
                        routine.poseidon_previous_state = routine.poseidon_end_state
                    except:
                        routine.poseidon_previous_state = routine.walk(Call)[0].poseidon_end_state
                        routine.poseidon_end_state = routine.poseidon_previous_state

            # gen summary
            meta = {
                'routines': {},
                'calls': {}
            }

            # fill routines
            cnt_streams = 0
            for name, routine in self.routines.items():
                entry = {
                    'prev': routine.poseidon_previous_state,
                    'end': routine.poseidon_end_state
                }
                meta['routines'][name] = entry
                cnt_streams += len(routine.poseidon_end_state['streams'])

            # fill calls
            for name, call in self.calls.items():
                entry = {
                    'prev': routine.poseidon_previous_state,
                    'end': routine.poseidon_end_state
                }
                meta['calls'][name] = entry

            cnt += 1
            if cnt > recurse:
                # final result
                self.meta = meta
                raise Exception("FAIL TO CONVERGE")
            
        # final result
        self.meta = meta
