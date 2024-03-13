import random
from typing import Tuple, Dict
import numpy as np
from functools import wraps


def remove_digits(func):
    """
    to represent viterbi decoding above the sequence - 1base=1state
    BS1S2S3 --> BSSS
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        res = func(*args, **kwargs)
        res = ''.join(i for i in res if not i.isdigit())
        return res
    return wrapper


def random_choice(prob_table: Dict[str, float]) -> str:
    """
    choose one random key weighted by values
    :param prob_table: {'hot': 0.7, 'cold': 0.2, 'warm': 0.1}
    :return: 'hot'
    """
    keys = list(prob_table.keys())
    values = list(prob_table.values())
    return random.choices(keys, weights=values, k=1)[0]


class State:
    def __init__(self, transitions: Dict[str, float], emissions=None):
        self.transitions = transitions
        self.emissions = emissions

    def generate_transition(self):
        return random_choice(self.transitions)

    def generate_emission(self):
        return random_choice(self.emissions)


class HMM:
    def __init__(self, states: Dict[str, State]):
        self.states = states
        self.fill_transitions()
        self.add_log_p()

    def fill_transitions(self):
        """
        If the transitions have not been transmitted, fills them with 0
        Example:
        model have states ['begin', 'background', 'site1', 'site2', 'site3', 'end']
        given State(transmissions={'background: 1}
        --> {'background': 1, 'begin': 0, 'site1': 0, ...}
        """
        all_states = list(self.states) + ['end']
        for state in self.states.values():
            for i in all_states:
                state.transitions.setdefault(i, 0)

    def add_log_p(self):
        """
        convert all probabilities to log
        in order not to get 0 when multiplying probabilities
        """
        for state in self.states.values():
            state.log_transitions = {i: np.log(j) for i, j in state.transitions.items()}
            if state.emissions:
                state.log_emissions = {i: np.log(j) for i, j in state.emissions.items()}


    def generate_sequence(self, length: int = None):
        """
        generate random sequence
        if length is given, generate sequence of that length
        if not, generate sequence of random length
        :return:
        """
        sequence = ""
        current_state = self.states['begin'].generate_transition()
        if length is None:
            while current_state != 'end':
                new_symbol = self.states[current_state].generate_emission()
                sequence += new_symbol
                current_state = self.states[current_state].generate_transition()
        else:
            for _ in range(length):
                new_symbol = self.states[current_state].generate_emission()
                sequence += new_symbol
                # transition
                new_state = self.states[current_state].generate_transition()
                # ignore end state, NOT THE MOST OPTIMAL WAY, BUT THE EASIEST FOR PROGRAMMING
                while new_state == "end":
                    new_state = self.states[current_state].generate_transition()
                current_state = new_state
        return sequence


    @remove_digits
    def viterbi(self, sequence: str) -> str:
        """
        performs viterbi decoding of given sequence
        """
        # initialize state_graph {'state1': [likelihoods], 'state2':[likelihoods] ...}
        states_graph = {key: [] for key in self.states.keys() if key not in ('begin', 'end')}

        # transitions from begin to the first state. u_first = a(beg, state) * e(state, first_symbol)
        for key in states_graph.keys():
            transition_b_key = self.states['begin'].log_transitions[key]
            emission_key_s0 = self.states[key].log_emissions.get(sequence[0], -np.inf)
            u = transition_b_key + emission_key_s0
            states_graph[key].append((u, 'begin'))

        # fill states graph using dynamic programing
        # go through each character of the string
        # `i` is the index of the previous state
        for i, symbol in enumerate(sequence[1:]):
            # choose the best path to each state and memorize previous state
            # u = max(u(i-1) * a(i-1, i) * e(state, symbol)
            for curr_state in states_graph.keys():
                u_to_choose = {}
                # consider each previous state
                for prev_state in states_graph.keys():
                    u_prev = states_graph[prev_state][i][0]
                    transition = self.states[prev_state].log_transitions[curr_state]
                    emission = self.states[curr_state].log_emissions.get(symbol, -np.inf)
                    u_to_choose[prev_state] = u_prev + transition + emission
                u_max = max(u_to_choose.values())
                state_from = max(u_to_choose, key=u_to_choose.get)
                states_graph[curr_state].append((u_max, state_from))

        # transitions to the end
        u_last_dict = {}
        for state in states_graph.keys():
            u_prev = states_graph[state][-1][0]
            transition = self.states[state].log_transitions['end']
            u_last_dict[state] = u_prev + transition
        u_last = max(u_last_dict.values())
        last_from = max(u_last_dict, key=u_last_dict.get)

        # backtracking to find optimal decoding
        i = -1
        decoding = ''
        prev_state = last_from
        while prev_state != 'begin':
            decoding += prev_state
            prev_state = states_graph[prev_state][i][1]
            i -= 1

        return decoding[::-1]


    def forward(self, sequence: str) -> Tuple[Dict[str, list], float]:
        """
        performs forward pass of forward-backward algorithm
        :return: state_graph with f values
        {'B': [0.01, 0.01, ...], 'S0': [0.0, 0.12, ...], ...}
        f_last - likelihood of the sequence - P(X)
        """
        # initialize state_graph {'state1': [likelihoods], 'state2':[likelihoods] ...}
        states_graph = {key: [] for key in self.states.keys() if key not in ('begin', 'end')}

        # transitions from begin to the first state. f_first = a(beg, state) * e(state, first_symbol)
        for key in states_graph.keys():
            transition_b_key = self.states['begin'].transitions[key]
            emission_key_s0 = self.states[key].emissions.get(sequence[0], 0)
            f = transition_b_key * emission_key_s0
            states_graph[key].append(f)

        # fill states graph using dynamic programing
        # go through each character of the string
        # `i` is the index of the previous state
        for i, symbol in enumerate(sequence[1:]):
            # for each current state sum f(i-1)*trans*emis
            for curr_state in states_graph.keys():
                f_to_sum = []
                # consider each previous state
                for prev_state in states_graph.keys():
                    f_prev = states_graph[prev_state][i]
                    transition = self.states[prev_state].transitions[curr_state]
                    emission = self.states[curr_state].emissions.get(symbol, 0)
                    f_to_sum.append(f_prev * transition * emission)
                f = sum(f_to_sum)
                states_graph[curr_state].append(f)

        # count f_last = likelihood of the sequence P(x)
        f_last_to_sum = []
        for state in states_graph.keys():
            f_prev = states_graph[state][-1]
            transition = self.states[state].transitions['end']
            f_last_to_sum.append(f_prev * transition)
        f_last = sum(f_last_to_sum)

        return states_graph, f_last

    def backward(self, sequence: str) -> Dict[str, list]:
        """
        performs backward pass of forward-backward algorithm
        :param sequence:
        :return:
        """
        # This dynamic programming goes from the end to the beginning, so it is different from viterbi and forward

        states_graph = {key: [0] * len(sequence) for key in self.states.keys() if key not in ('begin', 'end')}

        # b for the latest symbol. b_last = a(last, end)
        for key in states_graph.keys():
            transition_key_end = self.states[key].transitions['end']
            states_graph[key][-1] = transition_key_end

        # main part
        # from second last to the first
        for i in range(len(sequence) - 2, -1, -1):
            # b(i) = sum( b(i+1) * a(i, i+1) * e(i+1))
            for curr_state in states_graph.keys():
                b_to_sum = []
                for next_state in states_graph.keys():
                    b_next = states_graph[next_state][i + 1]
                    transition = self.states[curr_state].transitions[next_state]
                    emission = self.states[next_state].emissions.get(sequence[i + 1], 0)
                    b_to_sum.append(b_next * transition * emission)
                b = sum(b_to_sum)
                states_graph[curr_state][i] = b

#        # likelihood of the sequence. P(X)
#        b_first_to_sum = []
#        for first_state in states_graph.keys():
#            b_first = states_graph[first_state][0]
#            transition = self.states['begin'].transitions[first_state]
#            emission = self.states[first_state].emissions.get(sequence[0], 0)
#            b_first_to_sum.append(b_first * transition * emission)
#        b_first = sum(b_first_to_sum)

        return states_graph

    def forward_backward(self, sequence: str) -> Dict[str, list]:
        """
        performs aposterior decoding using forward-backward algorithm
        :return: Probabilites from FB for each state
        {'B': [0.01, 0.01, ...], 'S0': [0.0, 0.12, ...], ...}
        """
        forward, p_seq = self.forward(sequence)
        backward = self.backward(sequence)
        fb_probas = {}
        for state in forward.keys():
            fb_probas[state] = [(f * b) / p_seq for f, b in zip(forward[state], backward[state])]
        return fb_probas


def create_restriction_model(site, p_site, p_site_to_end=0.01, p_bg_to_end=0.01) -> HMM:
    """
    creates restriction site model with given params
    :param site: sequence of the restriction site
    :param p_site: probability of transition Background --> Site
    :param p_site_to_end: probability of transition Site --> End
    :param p_bg_to_end: probability of transition Background --> End
    """
    # create site_states
    site_states = {}
    for i, base in enumerate(site[:-1]):
        site_states[f"S{i}"] = State(transitions={f"S{i + 1}": 1},
                                     emissions={base: 1})
    site_states[f"S{len(site) - 1}"] = State(transitions={"end": p_site_to_end, "B": 1 - p_site_to_end},
                                             emissions={f"{site[-1]}": 1})
    # create begin and background
    states = {"begin": State(transitions={"B": 1}),
              "B": State(transitions={"S0": p_site, "end": p_bg_to_end, "B": 1 - p_site - p_bg_to_end},
                         emissions={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25})
              }

    states = {**states, **site_states}
    return HMM(states)

