#!/usr/bin/env python3
'''Generate modified Hamming code by dynamic programing.
The method is designed for small numbers of on-bits such as 4 on-bits.
Author: Xingjie Pan
'''

def H_dist(code1:set, code2:set):
    return len(code1.union(code2)) - len(code1.intersection(code2))  

def generate_modified_hamming_codes(length:int, n_on_bits:int, dist_cut:int):
    '''Generate a complete set of modified hamming codes.
    A complete set is a set that adding any extra code would
    break the required distance rule.

    Arguments:
        length: length of a code.
        n_on_bits: number of bits that are 1 for each code.
        dist_cut: the minimum distance between codes.
   
    Return:
        code_list: a list of the complete code set.
    '''
    assert(length > 0 and 0 <= n_on_bits <= length)

    code_list = []
    chosen_on_bits = set()

    gen_MHC_recursive(length, n_on_bits, dist_cut, code_list, chosen_on_bits)

    print(f'Generated a list of {len(code_list)} codes.')
    return code_list
   
def gen_MHC_recursive(length:int, n_on_bits:int, dist_cut:int, code_list:list, chosen_on_bits:set,
        print_depth:int=1):
    '''Generate modified hamming codes recursively.'''
    
    # If enough on-bits are chosen, add the code to the list
    if len(chosen_on_bits) == n_on_bits:
        code_list.append(chosen_on_bits)
        return

    # Find the candidates for the next on-bit

    on_bit_candidates = [i for i in range(length) if not i in chosen_on_bits]

    # The number of on bits after this round of addition

    n_rest_on_bits = n_on_bits - len(chosen_on_bits) - 1

    # Test all on-bit candidates

    for on_bit in on_bit_candidates:
        new_chosen_on_bits = chosen_on_bits.copy()
        new_chosen_on_bits.add(on_bit)

        if len(new_chosen_on_bits) <= print_depth:
            print(f'Testing partial set {new_chosen_on_bits}')

        ok_to_proceed = True

        # Compare the partly determined code to the existing codes

        for existing_code in code_list:
            Hd = H_dist(existing_code, new_chosen_on_bits)

            # Give up the partly determined code if the max possible distance is still smaller than the cutoff
            if Hd + n_rest_on_bits < dist_cut:
                ok_to_proceed = False
                break

        if ok_to_proceed:
            gen_MHC_recursive(length, n_on_bits, dist_cut, code_list, new_chosen_on_bits)

if __name__ == '__main__':
    code_list = generate_modified_hamming_codes(16, 4, 4)
    print(code_list)
