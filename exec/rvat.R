# collect args
collected_args <- rvat:::collect_args()

# run cli
rvat:::rvat_cli(
  args = collected_args[["args"]],
  args_raw = collected_args[["args_raw"]]
)