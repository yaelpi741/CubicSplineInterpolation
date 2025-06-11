def cubic_spline_interpolation(x_vals, y_vals, x_target):
    n = len(x_vals)

    # Input validation
    if n != len(y_vals):
        raise ValueError("The number of X values must match the number of Y values.")
    if n < 3:
        raise ValueError("At least 3 points are required for cubic spline interpolation.")
    if any(x_vals[i] >= x_vals[i + 1] for i in range(n - 1)):
        raise ValueError("X values must be sorted in strictly increasing order with no duplicates.")
    if x_target < x_vals[0] or x_target > x_vals[-1]:
        raise ValueError(f"Target point ({x_target}) is out of bounds: {x_vals[0]} to {x_vals[-1]}.")

    # Step 1: Calculate h and alpha
    h = [x_vals[i + 1] - x_vals[i] for i in range(n - 1)]
    alpha = [0] * n
    for i in range(1, n - 1):
        alpha[i] = (3 / h[i]) * (y_vals[i + 1] - y_vals[i]) - (3 / h[i - 1]) * (y_vals[i] - y_vals[i - 1])

    # Step 2: Solve tridiagonal system
    l = [1] + [0] * (n - 1)
    mu = [0] * n
    z = [0] * n
    for i in range(1, n - 1):
        l[i] = 2 * (x_vals[i + 1] - x_vals[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]

    l[-1] = 1
    z[-1] = 0

    # Step 3: Back-substitution
    c = [0] * n
    b = [0] * (n - 1)
    d = [0] * (n - 1)
    for j in range(n - 2, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]
        b[j] = (y_vals[j + 1] - y_vals[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3
        d[j] = (c[j + 1] - c[j]) / (3 * h[j])

    # Step 4: Evaluate the spline at the target point
    for i in range(n - 1):
        if x_vals[i] <= x_target <= x_vals[i + 1]:
            dx = x_target - x_vals[i]
            y_interp = y_vals[i] + b[i] * dx + c[i] * dx ** 2 + d[i] * dx ** 3
            return y_interp

    raise RuntimeError("No suitable interval found for x_target. Unexpected error.")


def main():
    try:
        print("Cubic Spline Interpolation")
        x_vals = [1, 2, 3, 4]       # X values
        y_vals = [1, 4, 9, 16]      # Y values
        x_target = float(input("Enter x value to interpolate: "))

        result = cubic_spline_interpolation(x_vals, y_vals, x_target)
        print(f"Estimated value at x = {x_target}: {result:.6f}")

    except ValueError as ve:
        print("ValueError:", str(ve))
    except RuntimeError as re:
        print("RuntimeError:", str(re))
    except Exception as e:
        print("Unexpected error:", str(e))

if __name__ == "__main__":
    main()
