import formStyles from '../../globalStyles/form.module.css'

const ReagentSelector = ({ ref }) => {
  return (
    <select className={formStyles.dropdown} ref={ref}>
      <option value=''>Select reagent ion</option>
      <option value=''>Any</option>
      <option value='NH4'>NH4+</option>
      <option value='NO'>NO+</option>
      <option value='H3O'>H3O+</option>
    </select>
  );
}

export default ReagentSelector;
