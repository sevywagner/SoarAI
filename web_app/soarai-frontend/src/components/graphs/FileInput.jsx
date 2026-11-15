import { useEffect, useState } from "react";
import styles from './file.module.css'
import ReagentSelector from '../styledComponents/ReagentSelector';

const FileInput = ({ reagentRef, onSubmit }) => {
  const [av, setAV] = useState();
  const [base, setBase] = useState();

  const checkFile = (file) => {
    if (!file) {
      return false;
    }
    
    if (file.name.split('.')[1] !== 'csv') {
      alert("File must be in csv format");
      return false;
    }
    return true
  }

  const avChangeHandler = (event) => {
    if (!checkFile(event.target.files[0])) {
      return;
    }
    setAV(event.target.files[0]);
  }

  const baseChangeHandler = (event) => {
    if (!checkFile(event.target.files[0])) {
      return;
    }
    setBase(event.target.files[0]);
  }

  const submitHandler = (e) => {
    e.preventDefault();

    if (!av || !base) {
      alert("Please updload both an mz_av and mz_base file");
      return;
    }

    onSubmit(av, base);
  }

  return (
    <div className={styles.wrapper}>
      <form className={styles.form} onSubmit={submitHandler}>
	<ReagentSelector ref={reagentRef} />
        <div className={styles.wrapper}>
          <div className={styles.input}>
            <label for="file1" className={av ? styles.green : styles.file}>MZ_AV</label>
            <input id="file1" type="file" name="mz_av" onChange={avChangeHandler} />
          </div>
          <div className={styles.input}>
            <label for="file2" className={base ? styles.green : styles.file}>MZ_Base</label>
            <input id="file2" type="file" name="mz_base" onChange={baseChangeHandler} />
          </div>
        </div>
        <button>Upload</button>
      </form>
    </div>
  );
}

export default FileInput;
