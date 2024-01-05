import numpy as np
import math
from scipy.optimize import linprog as linprog_scipy

class utils:

    def gauss_method(
            tableau: np.matrix[float],
            r: int,
            s: int
        ) -> np.matrix[float]:
            tableau[r] /= tableau[r,s]
            args = np.delete(np.arange(tableau.shape[0]), r)
            for arg in args:
                tableau[arg] -= tableau[arg][s] * tableau[r]
            return tableau

    def choose_s(tableau, case_condition = 1):
         z_line = tableau[0]
         def case_1():
              res = np.where(z_line < 0, z_line, np.inf).argmin()
              return res
         def case_2():
              s_candidates = np.where(z_line < 0)[0]
              return s_candidates[np.argmin(z_line[s_candidates])]
         def case_3():
              #q_i = tableau[1:,0]
              Tetta = np.zeros(tableau.shape[1])
              for j in range(1, tableau.shape[1]):
                   Tetta[j] = np.min(np.array([tableau[i,0] / tableau[i,j] for i in range(1, tableau.shape[0]) if tableau[i,j] > 0]))
              s = np.where(z_line < 0, Tetta, np.inf).argmax()
              return s
         def switch_case(case_number):
          if case_number == 1:
              return case_1()
          elif case_number == 2:
              return case_2()
          elif case_number == 3:
              return case_3()
          else:
              return
         return switch_case(case_condition)
    
    def cannonize(problem: dict):
        '''
        problem: 
            c: z_line
            A: matrix
            b: values
            type: max / min
            signs: (>=, =, <=)
        '''
        cn = np.copy(problem['c'])
        An = np.copy(problem['A'])
        bn = np.copy(problem['b'])
        signs = problem['signs']
        p_type = problem['type']
        basis = np.array([np.nan]*len(signs))

        if not all(s in ('<=', '=', '>=') for s in signs):
            raise Exception('Проверь знаки! У тебя ошибка')
        if problem['type'] not in ('max', 'min'):
            raise Exception('Тип задачи не соответствует симплекс-методу!')
        
        if p_type == 'min':
            cn = -cn
        
        aux_vars = []
        count =0
        for i in range(len(signs)):
            if signs[i] == '<=':
                aux_vars.append(1)
                count +=1
                basis[i] = An.shape[1] + i+1
            elif signs[i] == '>=':
                aux_vars.append(-1)
                count +=1
            else:
              aux_vars.append(0)
        aux_vars = np.array(aux_vars)
        if not np.all(aux_vars == 0):
          columns = np.array([np.zeros(count)]*len(signs))
          for j in range(count):
            for i in range(len(signs)):
              if aux_vars[i] == 0: pass
              else:
                columns[i,j] = aux_vars[i]
                aux_vars[i] = 0
                j +=1
          An = np.hstack((An,columns))
          c_vars = np.zeros(aux_vars.shape[0])
          cn = np.append(cn, c_vars)
        b = np.insert(bn,0,0)
        while np.isnan(basis).any():
          None_indices = np.where(np.isnan(basis))[0]
          for i in None_indices:
            column = np.zeros(len(signs))
            column[i] = 1

            An = np.hstack((An, column.reshape(-1, 1)))
            basis[i] = An.shape[1]
        c_vars = np.zeros(An.shape[1] - problem['A'].shape[1])
        cn = np.append(cn, c_vars)
        return cn, An,b, basis



class _SymplexTable:

    def __init__(
            self,
            problem: dict
        ):
            self.table = NotImplemented
            self._has_optimal_solution = False
            self._was_solved = False
            self.basis = NotImplemented
            z_line, A, b, basis = utils.cannonize(problem)
            self._set_table(z_line, A, b, basis)
    
    def _set_table(self,func, A_matrix, B, basis):
        '''
        :параметр: np.ndarray func: Коэффициенты уравнения,
        :параметр: np.ndarray A_matrix: Левая часть системы ограничений,
        :параметр: np.ndarray B: Правая часть системы ограничений.
        :возврат: None
        '''

        self._check_input_data(func, A_matrix, B,basis)
        self.table = np.zeros([size+1 for size in A_matrix.shape], dtype= np.float64)
        # Канонизация введенных параметров

        self.table[1:, 1:] = A_matrix
        self.table[:, 0] = B
        self.table[0,1:] = func
        self.basis = basis
        print(self.table)
        print(f"start_basis: {basis}")

    @staticmethod
    def _check_input_data(
            func: np.ndarray,
            a_matrix: np.ndarray,
            b_vec: np.ndarray,
            basis: np.ndarray 
    ) -> None:
        pass 
    

    def iterate(self)->dict:

        # Критерий выхода из алгоритма
        def criterion(z_line):
            return np.any(z_line < 0)
        
        q = self.table
        iter = 0
        while criterion(q[0,1:]):
            
            #Выбор s: q[0,s] < 0:
            s = utils.choose_s(q)

            if not np.any(q[:,s] >0):
                # Если все элементы в столбце s неположительные, то значение целевой функции неограничено...
                self._has_optimal_solution = False
                break

            else:
                q_0 = q[1:,0]
                q_s = q[1:,s]               
                relations = np.divide(q_0, q_s, out = np.full_like(q_s, np.inf), where= (q_s != 0))
                r = np.argmin(relations)
                q = utils.gauss_method(q,r,s)
                self.basis[r-1] = s
            iter +=1
            if iter > 50:
              print(iter)
              break
        return {'table': self.table,
                'basis': self.basis,
                'Solution': self._has_optimal_solution,
                'solved': self._was_solved}