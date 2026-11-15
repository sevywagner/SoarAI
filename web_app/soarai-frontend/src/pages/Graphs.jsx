import FileInput from "../components/graphs/FileInput";
import { useEffect, useState, useRef } from "react";
import formStyles from './../globalStyles/form.module.css'
import wrapperStyles from './../globalStyles/wrappers.module.css'
import textStyles from '../globalStyles/text.module.css'
import styles from './css/graphs.module.css'
import BackButton from "../components/styledComponents/BackButton";

const Graphs = () => {
  const [img, setImg] = useState();
  const [time, setTime] = useState();
  const [inputMode, setInputMode] = useState(true);
  const [loading, setLoading] = useState(false);
  const [preds, setPreds] = useState();
  const [errorMsg, setErrorMsg] = useState();
  const reagentRef = useRef(null);
  
  const submitHandler = (av, base) => {
    setErrorMsg(null);
    const form = new FormData();

    if (!av || !base) {
      setErrorMsg("Upload CSVs for both av and base");
      return;
    }
    
    form.append('av', av);
    form.append('base', base);
    form.append('reagentIon', reagentRef.current.value);

    let hasErr = false;
    let errStatus;
    setLoading(true);
    fetch('http://localhost:18080/process-graph', {
      method: 'POST',
      body: form
    }).then((res) => {
      if (!res.ok) {
	hasErr = true;
	errStatus = res.status;
      }
      
      return res.text();
    }).then((data) => {
      if (hasErr) {
	setErrorMsg(data);
	console.log(status + ": " + data);
	setLoading(false);
	return;
      }

      setLoading(false);
      
      setPreds(URL.createObjectURL(new Blob([data], { type: 'text/plain' })));
      setInputMode(false);
    }).catch((err) => {
      setLoading(false);
      console.log(err);
      setErrorMsg("Error occurred. Please try again.");	
    });
  }

  const toggleInputMode = () => {
    setInputMode((prevState) => !prevState);
  }

  useEffect(() => {
    if (inputMode) {
      setPreds(null);
    }
    
    setErrorMsg(null);
  }, [inputMode, preds])

  return (
    <div>
      <div className={wrapperStyles.three}>
        <div className={wrapperStyles['three-small']}>
          {!inputMode && img && <BackButton onClick={toggleInputMode}>Back</BackButton>}
        </div>
        <div className={wrapperStyles['three-big']}>
          <p className={textStyles.title}>Upload Graph</p>
        </div>
        <div className={wrapperStyles['three-small']}></div>
      </div>
      {errorMsg && <p className={formStyles.error}>{errorMsg}</p>}
      {inputMode && <FileInput reagentRef={reagentRef} onSubmit={submitHandler} />}
      {!inputMode &&
       <div className={wrapperStyles.center}>
         {!inputMode && preds &&
	  <a href={preds} className={styles.option} download="output.txt">Download Results</a>}
       </div>}
      {loading && <div className={wrapperStyles.center}> <p>Loading...</p></div>}
      {!inputMode && img && <img src={`data:image/png;base64,${img}`} />}
    </div>
  );
}

export default Graphs;
