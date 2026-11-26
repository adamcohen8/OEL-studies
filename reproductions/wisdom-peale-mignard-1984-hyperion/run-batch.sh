set -eu
which_run=$1
case "$which_run" in
 production)
 set -- production --out outputs/production-c16
 for theta in 0 0.7853981633974483; do
   for w in -0.5 0 0.5 1 1.5 2 2.65; do set -- "$@" --seed "$theta" "$w"; done
 done
 set -- "$@" --orbits 64 --steps-per-orbit 128 --orbit-integrator rkf78 --event-every 4 --no-lyap --no-periodic
 ;;
 refine) set -- refine --reference outputs/production-c16 --out outputs/refine-c16 --orbits 64 --steps-per-orbit 256 --orbit-integrator rkf78 --event-every 4 --no-lyap --no-periodic ;;
 targeted) set -- production --out outputs/targeted-c16 --seed 0 0 --seed 0 0.5 --seed 0 1 --seed 0 1.5 --seed 0 2.65 --seed 0.7853981633974483 2.65 --orbits 96 --steps-per-orbit 256 --orbit-integrator rkf78 --event-every 1 --lyap --periodic --search-min -1.15 --search-max 2.9 --search-grid 33 ;;
 long) set -- production --out outputs/long-c16 --seed 0 1.5 --orbits 384 --steps-per-orbit 256 --orbit-integrator rkf78 --event-every 1 --lyap --no-periodic ;;
 halfstep) set -- production --out outputs/halfstep-c16 --seed 0 1.5 --orbits 192 --steps-per-orbit 512 --orbit-integrator rkf78 --event-every 1 --lyap --periodic --search-min -1.15 --search-max 2.9 --search-grid 33 ;;
 *) exit 2 ;;
esac
"${STUDY_PYTHON:-python}" replay.py "$@"
