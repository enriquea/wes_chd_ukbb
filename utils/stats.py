# eam
# 2021-02-03

import hail as hl


# run two-tailed fisher exact test
def compute_fisher_exact(tb: hl.Table,
                         n_cases_col: str,
                         n_control_col: str,
                         total_cases_col: str,
                         total_controls_col: str,
                         correct_total_counts: bool,
                         extra_fields: dict,
                         root_col_name: str = 'fet') -> hl.Table:
    """
    Perform two-sided Fisher Exact test. Add extra annotations (if any)

    :param tb: Hail Table
    :param n_cases_col: field name with number of (affected) cases
    :param n_control_col: field name with number of (affected) control
    :param total_cases_col: field name with total number of cases
    :param total_controls_col: field name with total number of controls
    :param correct_total_counts: should the total numbers (case/control) be corrected to avoid duplicated counting?
    :param root_col_name: field to be annotated with test results
    :param extra_fields: Extra filed (must be a dict) to be annotated
    :return: Hail Table with Fisher Exact test results.
    """
    # compute fisher exact
    if correct_total_counts:
        fet = hl.fisher_exact_test(c1=hl.int32(tb[n_cases_col]),
                                   c2=hl.int32(tb[n_control_col]),
                                   c3=hl.int32(tb[total_cases_col]) - hl.int32(tb[n_cases_col]),
                                   c4=hl.int32(tb[total_controls_col]) - hl.int32(tb[n_control_col]))
    else:
        fet = hl.fisher_exact_test(c1=hl.int32(tb[n_cases_col]),
                                   c2=hl.int32(tb[n_control_col]),
                                   c3=hl.int32(tb[total_cases_col]),
                                   c4=hl.int32(tb[total_controls_col]))

    tb = (tb
          .annotate(**{root_col_name: fet})
          .flatten()
          )

    if len(extra_fields) == 0:
        return tb
    else:
        return tb.annotate(**extra_fields)


# run MT rows logistic regression
def logistic_regression(mt: hl.MatrixTable,
                        x_expr: str,
                        response: str,
                        covs: list,
                        pass_through: list,
                        extra_fields: dict,
                        add_odd_stats: bool = True) -> hl.Table:
    """
    Perform a logistic-regression test (use by default Wald test).

    :param mt: Hail MatrixTable
    :param x_expr: the genotype field name (numeric expression)
    :param response: binary response
    :param covs: list of covariates to be included in the test
    :param pass_through: list of extra fields to keep in the output
    :param extra_fields: extra field to annotated (expected a dict)
    :param add_odd_stats: compute odds from logistic regression stats
    :return: Hail Table with logistic regression test results
    """
    # parsing covariates list
    if len(covs) >= 1:
        covs = [1] + [mt[field] for field in covs]
    else:
        covs = [1]

    tb_stats = hl.logistic_regression_rows(y=mt[response],
                                           x=mt[x_expr],
                                           covariates=covs,
                                           pass_through=pass_through,
                                           test='wald')

    if add_odd_stats:
        # Compute Odds ratio and 95% confidence interval from logistics regression stats
        tb_stats = tb_stats.annotate(odds_ratio=hl.exp(tb_stats.beta),
                                     lower_ci_95=hl.exp(tb_stats.beta - 1.96 * tb_stats.standard_error),
                                     upper_ci_95=hl.exp(tb_stats.beta + 1.96 * tb_stats.standard_error))

    # add column with additional information
    if len(extra_fields) == 0:
        return tb_stats
    else:
        return tb_stats.annotate(**extra_fields)
