"""Independent dense partial-pivot elimination, no numerical dependencies."""
import argparse
import json
import math
import subprocess


def reference(nx, ny):
    hx, hy = 1.5 / nx, 0.75 / ny
    n = (nx - 1) * (ny - 1)
    a = [[0.0] * n for _ in range(n)]
    rhs = [0.0] * n
    def k(i, j):
        return 1.0 + i * hx + 0.5 * j * hy
    def boundary(i, j):
        x, y = i * hx, j * hy
        return x + 2 * y + 0.1 * x * y
    for j in range(1, ny):
        for i in range(1, nx):
            row = (j - 1) * (nx - 1) + i - 1
            x, y = i * hx, j * hy
            rhs[row] = 1 + x * y + x * x
            a[row][row] = 0.3
            for ii, jj, h in [(i-1,j,hx),(i+1,j,hx),(i,j-1,hy),(i,j+1,hy)]:
                weight = (2 * k(i,j) * k(ii,jj) / (k(i,j)+k(ii,jj))) / h**2
                a[row][row] += weight
                if 0 < ii < nx and 0 < jj < ny:
                    a[row][(jj-1)*(nx-1)+ii-1] -= weight
                else:
                    rhs[row] += weight * boundary(ii,jj)
    # Independently constructed dense operator; do not use exported MoonBit COO.
    for col in range(n):
        pivot = max(range(col,n), key=lambda r: abs(a[r][col]))
        a[col], a[pivot] = a[pivot], a[col]
        rhs[col], rhs[pivot] = rhs[pivot], rhs[col]
        assert abs(a[col][col]) > 1e-14
        for row in range(col+1,n):
            factor = a[row][col] / a[col][col]
            for c in range(col+1,n):
                a[row][c] -= factor * a[col][c]
            rhs[row] -= factor * rhs[col]
    u = [0.0] * n
    for row in reversed(range(n)):
        u[row] = (rhs[row] - sum(a[row][c]*u[c] for c in range(row+1,n))) / a[row][row]
    full = []
    for j in range(ny+1):
        for i in range(nx+1):
            full.append(boundary(i,j) if i in (0,nx) or j in (0,ny) else u[(j-1)*(nx-1)+i-1])
    return full


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--target', default='js')
    args = parser.parse_args()
    run = subprocess.run(['moon','run','cmd/reference','--target',args.target], check=True, stdout=subprocess.PIPE, universal_newlines=True)
    samples = [json.loads(s) for s in run.stdout.splitlines() if s.startswith('{')]
    assert len(samples) == 12, 'missing reference samples'
    worst = 0.0
    for sample in samples:
        expected = reference(sample['nx'],sample['ny'])
        assert len(expected) == len(sample['solution'])
        for actual, wanted in zip(sample['solution'],expected):
            assert math.isfinite(actual)
            error = abs(actual-wanted)
            worst = max(worst,error)
            assert error < 1e-9, (sample['nx'],sample['ny'],error)
    print('differential: 12/12 grids passed; max absolute nodal error %.3g' % worst)


if __name__ == '__main__':
    main()
